template<typename T>
T pbc(T input, T size)
    {
    if (input < 0)
        input += size;
    if (input >= size)
        input -= size;
    return input;
    }
    
    
template<typename T>
T ipbc(T pos, T box_half, T com) // normally pos[dim], box_half[dim], center of mass [dim]
{   
    if (com < box_half)
    {
        if ((pos-com)*(pos-com) > box_half*box_half)
        {   
            pos -= 2*box_half;
        }
    }
    
    if (com > box_half)
    {
        if ((pos-com)*(pos-com) > box_half*box_half)
        {
            pos += 2*box_half;
        }
    }
    return pos;
}
