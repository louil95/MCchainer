class Trig
{
private:
    // for rotation
    mutable double cosangle;
    mutable double min_cosa;
    mutable double sinangle;
  
public:
    void prepare_rodriguez(double angle)
    {
        cosangle = cos(angle);
        min_cosa = 1 - cosangle;
        sinangle = sin(angle);
    }
    
    void get_cross(double cross[3], double vec_0[3], double vec_1[3])
    {
        cross[0]=vec_0[1]*vec_1[2] - vec_0[2]*vec_1[1];
        cross[1]=vec_0[2]*vec_1[0] - vec_0[0]*vec_1[2];
        cross[2]=vec_0[0]*vec_1[1] - vec_0[1]*vec_1[0];    
    }
    double periodic_dist(double dist, int dim)
    {
        if (dist > box_size[dim]/2)
            dist -= box_size[dim];
        else if ((dist < -box_size[dim]/2))
            dist += box_size[dim];
        return dist;
    }
        
    double get_dist(double pos_0[3], double pos_1[3])
    {
        double dist;
        double sq_dist{0.};
        
        for (int dim = 0; dim < 3; ++dim)
            {
            dist = pos_0[dim] - pos_1[dim];
            dist = periodic_dist(dist, dim);
            sq_dist += dist*dist;
            }
        return sq_dist;
    }
    double get_angle(double a[3], double b[3], double c[3], double com[3])
    {   
        double ca[3];
        double cb[3];
        double angle;
        for (int dim = 0; dim<3; ++dim)
        {
            ca[dim] = ipbc(a[dim], box_size[dim]/2 ,com[dim]) - ipbc(c[dim], box_size[dim]/2,com[dim]);
            cb[dim] = ipbc(b[dim], box_size[dim]/2 ,com[dim]) - ipbc(c[dim], box_size[dim]/2,com[dim]);
        }
        normalize(ca);
        normalize(cb);
        angle = dot_prod(ca,cb);
        angle = acos(angle);
        
        return angle;
    }
    
    double get_tors(double r_0[3], double r_1[3], double r_2[3], double r_3[3], double com[3])
    {   
        double tors{0};
        double cross_0[3];
        double cross_1[3];
        double r_10[3];
        double r_21[3];
        double r_32[3];
        
        for (int dim = 0; dim<3; ++dim)
        {
            r_10[dim] = ipbc(r_1[dim], box_size[dim]/2 ,com[dim]) - ipbc(r_0[dim], box_size[dim]/2,com[dim]);
            r_21[dim] = ipbc(r_2[dim], box_size[dim]/2 ,com[dim]) - ipbc(r_1[dim], box_size[dim]/2,com[dim]);
            r_32[dim] = ipbc(r_3[dim], box_size[dim]/2 ,com[dim]) - ipbc(r_2[dim], box_size[dim]/2,com[dim]);
        }
        get_cross(cross_0, r_32, r_21);
        get_cross(cross_1, r_21, r_10);
        
        normalize(cross_0);
        normalize(cross_1);
        
        tors = dot_prod(cross_0,cross_1);
        tors*= (-1);
        return tors;       
    }
    
    void normalize(double vec[3])
    {   
        double norm{0.};
        
        for (int dim = 0; dim<3; ++dim)
        {
            norm += vec[dim] * vec[dim];
        }
        
        norm = sqrt(norm);
        
        for (int dim = 0; dim<3; ++dim)
        {
            vec[dim] /= norm;
        }
    }
    
    
    double dot_prod(double pos_0[3], double pos_1[3])
    {   
        double dot_product{0.};
        for (int dim=0; dim<3; ++dim)
            {
                dot_product += pos_0[dim] * pos_1[dim];
            }
        return dot_product;
    }
    
    void rodriguez(double k[3], double a[3])
    {   
        normalize(k);
        double k_dot_a = dot_prod(a,k);
        double a_[3]{a[0], a[1],a[2]};
        a[0] = a_[0] * cosangle + (k[1]*a_[2] - k[2]*a_[1]) * sinangle + k[0] * k_dot_a * min_cosa;    
        a[1] = a_[1] * cosangle + (k[2]*a_[0] - k[0]*a_[2]) * sinangle + k[1] * k_dot_a * min_cosa;   
        a[2] = a_[2] * cosangle + (k[0]*a_[1] - k[1]*a_[0]) * sinangle + k[2] * k_dot_a * min_cosa;
    
    }
    
    void center_of_mass(double center_of_mass[n_molecules][3], Segment molecules[n_molecules][n_segments])
    {   
        double theta_i;
        double cos_average;
        double sin_average;
        double factor[3];
        
        factor[0] = 2 * M_PI / box_size[0];
        factor[1] = 2 * M_PI / box_size[1];
        factor[2] = 2 * M_PI / box_size[2];

        for (int dim=0; dim<3; ++dim)
        {
            for (int mol=0; mol<n_molecules; ++mol)
            {   
                cos_average = 0.;
                sin_average = 0.;
                for (int seg=0; seg<n_segments; ++seg)                   
                {   
                    theta_i = molecules[mol][seg].C[dim];
                    theta_i *= factor[dim];
                    
                    cos_average += cos(theta_i);
                    sin_average += sin(theta_i);                   
                }
                cos_average /= n_segments;
                sin_average /= n_segments;
                
                center_of_mass[mol][dim] = atan2(-sin_average, -cos_average) + M_PI;
                center_of_mass[mol][dim] /= factor[dim];
            }
        }        
    }
    
    void one_com(double com[3], Segment molecules[n_molecules][n_segments], int mol)
    {   
        double theta_i;
        double cos_average;
        double sin_average;
        double factor[3];
        
        factor[0] = 2 * M_PI / box_size[0];
        factor[1] = 2 * M_PI / box_size[1];
        factor[2] = 2 * M_PI / box_size[2];

        for (int dim=0; dim<3; ++dim)
        {
            cos_average = 0.;
            sin_average = 0.;
            for (int seg=0; seg<n_segments; ++seg)                   
            {   
                theta_i = molecules[mol][seg].C[dim];
                theta_i *= factor[dim];
                
                cos_average += cos(theta_i);
                sin_average += sin(theta_i);                   
            }
            cos_average /= n_segments;
            sin_average /= n_segments;
                               
            com[dim] = atan2(-sin_average, -cos_average) + M_PI;
            com[dim] /= factor[dim];
        }
              
    }
    

    
    
};