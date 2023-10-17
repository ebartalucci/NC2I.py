# To Do List:

## Bugs to fix in scan_distances.py
* ~~at the moment it displaces the structure laterally and not along the centroid axis. i think we have a problem in the iterations, it might do some weird update in the displacement_direction variable or is it a fragments_displacer problem?~~
* displacement now happens along the centroid vector but it doesnt work all the time, sometimes gives crap, why??
* Error handling and log files: include folder - they need to update for each iteration
* ~~displacement steps increases for each iteration~~
* ~~code pauses run at each generated figure~~
* ~~check if coordinates are correct~~
* ~~correct these FUCKING relative paths~~
* ~~save output structures in correct subfolder~~
* ~~correct all paths to relative paths~~
* ~~at the end each displaced structure contains only the displaced fragment without the fixed one~~
* ~~atoms names are reshuffled and called Xx after displacement, keep them fixed to their original identity~~
* the code is not scalable for arbitrary number of fragments, can be a problem for multiple structures
* ~~need to merge the coordinates with the respective atom identities. at the moment it just give a X~~
* plot atom identities on the topology - the variable is atom_identities = atom_data[0]
* define a general fragment color but also would be cool to have atom based colors in the plot, for ease of representation
* keep in mind that the displacement direction will be fucked up if one introduces more than two fragments, since it cant displace them along the line connecting them anymore -> fixable by displacing fragment coordinates rather then centroids
* calculate and write to displaced file the distance between the centroid of the fixed fragment and each atom of the displaced fragment, this is needed for plotting the shifts vs distance data


## Structure
* ~~maybe make it modular so that we can reuse some functions into other implementations~~
* version control

## Features to implement
* dissociation limit as threshold for the foor loop - do i need to calculate it numerically or can i just extract it from the DFT optimization output file? this is apparently not so trivial
* write python output file with all infos of the run
* allwo displacement direction to go back some angstrom to access distances within the van der Waals sphere, i.e. shorter than the distance computed from the geometry optimization
* ~~number of K-means clusters = number of centroids to have = number of actual molecular fragments~~

## Style and appearance
* Dynamically updating figure with subpanels displaying:
  - initial topology
  - initial + displaced coordinates, with shaded structures in the middle
  - output summary including number and names of generated structures, date, running summary and description of input output files, all saved in a log file that is updated every run

## Long term
* can i machine learn the shift like this?
* include a MySQL database connection for future storage and machine learning of the shifts
* update so that also the output of the DFT calculations is displaced on the figure, after computation
* .exe file so that installation is easy for all

