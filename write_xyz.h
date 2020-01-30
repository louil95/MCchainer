#include <iostream>
#include <fstream>

void write_xyz(Segment (molecules[n_molecules][n_segments]), double com[n_molecules][3],std::string traj)
{   
    std::ofstream xyz;
  	xyz.open (traj, std::fstream::app);

  	    //xyz << "\n";
  	    xyz << n_molecules * (n_segments) << "\n" << "\t" << "\n";
  	    xyz.precision(4);
  	    
  	    for (int mol=0; mol < n_molecules; ++mol)
  	    {
  	        for (int seg=0; seg < n_segments; ++seg)
  	        {   
  	            xyz << "C "<< ipbc(molecules[mol][seg].C[0], box_size[0]/2, com[0][0]) - com[0][0]
  	                << "\t" << ipbc(molecules[mol][seg].C[1], box_size[1]/2., com[0][1]) - com[0][1] << "\t"
  	                << ipbc(molecules[mol][seg].C[2], box_size[2]/2, com[0][2]) - com[0][2] << "\n";
  	        }
  	    }
  	
	
	xyz.close();
}