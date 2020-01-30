#include <string>
#include <iostream>
#include <fstream>

std::string append_string(int count_atom, int seg,std::string atom_name, double x, double y, double z)
{   
    int digit;
    float remainder;
    int subst;
    int int_pos;
    
    std::string line{"ATOM  "};
    
    int string_size = ceil(float(count_atom)/100.);
    for (digit=0; digit < 5 - string_size; ++digit)
        line.append(" ");

    line.append(std::to_string(count_atom));
    atom_name.append(std::to_string(seg+1));
    string_size = atom_name.size();
    for (digit=0; digit < 4 - string_size; ++digit)
        line.append(" ");    
    line.append(atom_name);
    line.append("  ");
    line.append("ETTL    1      ");
    
    double xyz_arr[3]{x,y,z};
    
    for (int dim = 0; dim<3; ++dim)
    {
        string_size = ceil(xyz_arr[dim]/10.);
        string_size += 4;
        //std::cout<<8-string_size<<"\n";
        for (digit=0; digit < 8 - string_size; ++digit)
            line.append(" ");
        
        int_pos = int(floor(xyz_arr[dim]));
        line.append(std::to_string(int_pos));
        line.append(".");
        remainder = xyz_arr[dim] - int_pos;

        subst = int(floor(remainder*10));
        line.append(std::to_string(subst));
        remainder -= float(subst)/10.;

        subst = int(floor(remainder*100));
        line.append(std::to_string(int(floor(subst ))));
        remainder -= float(subst)/100.;

        subst = int(floor(remainder*1000));
        line.append(std::to_string(int(floor(remainder*1000))));
    }
    line.append("  1.00  0.00              \n");
    return line;   
}

void write_pdb(Segment (molecules[n_molecules][n_segments]))
{
    std::ofstream pdb;
  	pdb.open ("traj.pdb", std::fstream::app);
  	std::string line;
  	    
  	    for (int mol=0; mol < n_molecules; ++mol)
  	    {    
  	        int count_atom{0};
  	        for (int seg=0; seg < n_segments ; ++seg)
  	        {   

  	            count_atom++;
  	            line = append_string(count_atom, seg,"C", molecules[mol][seg].C[0], molecules[mol][seg].C[1] , molecules[mol][seg].C[2]);
  	            //pdb << "ATOM      " << count_atom << " " << "C" << seg+1 << "\tBUTL " << "1" << " " << molecules[mol][seg].C[0] << " " << molecules[mol][seg].C[1] << " " << molecules[mol][seg].C[2] << " 1.00 0.00\n";
  	            pdb << line;
  	        }
  	        for (int seg=0; seg < n_segments ; ++seg)
  	        {
  	            count_atom++;
  	            line = append_string(count_atom, seg,"1H", molecules[mol][seg].H0[0], molecules[mol][seg].H0[1] , molecules[mol][seg].H0[2]);
  	            //pdb << "ATOM      " << count_atom << " " << "1H" << seg+1 <<"\tBUTL " << "1" << " " << molecules[mol][seg].H0[0] <<" " << molecules[mol][seg].H0[1] << " " << molecules[mol][seg].H0[2] << "\t1.00\t0.00\n";
  	            pdb << line;
  	            count_atom++;
  	            line = append_string(count_atom, seg,"2H", molecules[mol][seg].H1[0], molecules[mol][seg].H1[1] , molecules[mol][seg].H1[2]);
  	            pdb << line;
  	            //count_atom++;
  	            //
  	            //pdb << "ATOM      " << count_atom << " " << "2H" << seg+1 <<"\tBUTL " << "1" << " " << molecules[mol][seg].H1[0] <<" " << molecules[mol][seg].H1[1] << " " << molecules[mol][seg].H1[2] << "\t1.00\t0.00\n";
  	            
  	            if (seg == 0 || seg == n_segments -1)
  	                count_atom++;
  	                line = append_string(count_atom, seg,"3H", molecules[mol][seg].H2[0], molecules[mol][seg].H2[1] , molecules[mol][seg].H2[2]);
  	                pdb << line;
  	                //pdb << "ATOM      " << count_atom << " " << "3H" << seg+1 <<"\tBUTL " << "1" << " " << molecules[mol][seg].H2[0] <<" " << molecules[mol][seg].H1[1] << " " << molecules[mol][seg].H1[2] << "\t1.00\t0.00\n";
  	        }
  	        pdb << "TER\n";
  	    }
  	    pdb << "ENDMDL\n\n\n";
  	       
  	          	
	pdb.close();
}