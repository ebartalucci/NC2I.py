################################################################
#          GUI for (NC)2Ipy software for noncov toolbox        #
#                                                              #
#            Ettore Bartalucci, v 0.1, 12.10.23 Aachen         #
################################################################

# NC2Ipy - NONCOV TOOLBOX

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QMenuBar, QMenu, QFileDialog, QTextEdit, QVBoxLayout, QWidget, QFrame, QSizePolicy, QHBoxLayout
from PyQt5.QtGui import QTextCursor
from PyQt5.QtCore import Qt
import importlib
import datetime

class TerminalWidget(QTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(False)
        self.setStyleSheet("background-color: black; color: white;")
        self.prompt = ">> "
        self.cursor = self.textCursor()
        self.insertPrompt()

        # Terminal logs
        self.log_file = open("gui/terminal_log.txt", "a")

        # Write the headline with date, time, and creator's name
        headline = f"-------- (NC)2I.py TERMINAL LOGS: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} created by Ettore Bartalucci -------- \n"
        self.log_file.write(headline)
        self.log_file.flush()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Return:
            self.insertPlainText("\n")
            self.processCommand()
            self.insertPrompt()
        else:
            super().keyPressEvent(event)

    def insertPrompt(self):
        self.insertPlainText(self.prompt)
        self.cursor.movePosition(QTextCursor.End)
        self.setTextCursor(self.cursor)

    def processCommand(self):
        # Get the command entered by the user
        command = self.toPlainText().split('\n')[-2].replace(self.prompt, '')

        # Clear terminal
        if command == 'clear':
            self.clear()
            self.insertPrompt()
        else:
            # Here, you can execute the command or perform any other desired action.
            # Write the terminal to log file
            self.log_file.write(f">> {command}\n")
            self.log_file.flush()

            # Display the entered command
            self.insertPlainText(f">> You entered: {command}\n")
            self.insertPrompt()

class NC2IpyGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("NC2Ipy - Scientific Software")
        self.setGeometry(100, 100, 800, 400)

        # Create a central widget
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        # Create a layout for the central widget
        main_layout = QHBoxLayout()

        # Create a menu bar
        menubar = self.menuBar()

        # File Menu
        file_menu = menubar.addMenu("File")
        open_action = QAction("Open", self)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        # Create additional menus and features
        self.create_menu("Structure search", ["PDBparser", "CCDCparser", "vdW_sphere_search", "HistoNCI"])
        self.create_menu("NONCOV", ["Geometrical analysis", "DFT on NMR", "Autoref shifts", "Tensorview", "RMSD NMR-DFT"])
        self.create_menu("NMR", ["Secondary structure CS", "Dipole-dipole calc", "Johnson-Bovey calc", "3D Clustering shifts", "PCA", "DIPSHIFT", "sedNMR"])
        self.create_menu("DFT", ["Generate structures", "CREST", "Tensors"])
        self.create_menu("Analysis", ["NCIplot", "Relaxation", "Orbitals"])
        self.create_menu("Tools", ["Generate Orca Input", "Generate CREST input"])

        # Create a frame for the terminal display
        terminal_frame = QFrame()
        terminal_frame.setFrameStyle(QFrame.Panel)
        terminal_frame.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)

        # Create the terminal widget
        self.terminal = TerminalWidget()

        # Add the terminal widget to the frame
        terminal_frame_layout = QVBoxLayout()
        terminal_frame_layout.addWidget(self.terminal)
        terminal_frame.setLayout(terminal_frame_layout)

        # Add the terminal frame to the main layout
        main_layout.addWidget(terminal_frame)

        central_widget.setLayout(main_layout)

    def create_menu(self, menu_name, features):
        new_menu = self.menuBar().addMenu(menu_name)
        for feature in features:
            action = QAction(feature, self)
            action.triggered.connect(self.call_external_function)
            new_menu.addAction(action)

    def call_external_function(self):
        # Get the text of the triggered action (feature)
        feature = self.sender().text()
        # Dynamically import the external function module
        try:
            module = importlib.import_module(feature)
            # Call the external function in the module (assumes the function has a predefined name)
            module.external_function()
        except ImportError as e:
            error_message = f"Module {feature}.py not found or does not contain the expected function."
            self.log_error(error_message)
        except Exception as e:
            error_message = f"Error in {feature}: {str(e)}"
            self.log_error(error_message)

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_name, _ = QFileDialog.getOpenFileName(self, "Open .xyz File", "", "XYZ Files (*.xyz)", options=options)
        # Replace this with the logic to handle the selected file

    def log_error(self, error_message):
        self.terminal.insertPlainText(error_message + "\n")

def main():
    app = QApplication(sys.argv)
    window = NC2IpyGUI()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
