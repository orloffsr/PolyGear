from math import cos, floor, ceil, pi



    
def ratio_search(teeth_min = 11, ring_max = 50,  
                planets_min = 3, planets_max = 12, 
                m2_min = 0.5, m2_max = 2.0,
                ratio_min = 1200.0, ratio_max = 1E10,
                results_max = 10000):
        """Performs exhaustive search and returns a list of tuples.
        each tuple consists of (Np, s1, p1, r1, s2, p2, r2, ratio, m2) 
        meeting the criteria given as input. 
        
        teeth_min = minimum number of teeth in any gear
        ring_max  = maximum number of teeth in ring gear
        planets_min and planets_max = min and max number of planets
        m2_min and m2_max = minimum and maximum module of second stage. m of first stage assumed to be 1.0.
        ratio_min and ratio_max = range of gear ratios desired.
        results_max = max number of results to return
        """
        results = []
        
        p_max = int((ring_max - teeth_min)/2)
        
        for Np in range(planets_min,planets_max+1):
            f = 0.5 if Np%2 ==0 else 1.0
            theta = pi*(1.0-2.0/Np)/2.0
            def zmin_fun(p):
                z = ceil(((p+1.0)*(1/cos(theta)-1)+p)/Np/f)
                return max(z,1)
                
            for p1 in range(teeth_min, p_max):
                zmin =  zmin_fun(p1) #planets can't intersect
                #for z1 in range(0,zmax):
                s1 = 0
                z1=-1
                while(s1+2*teeth_min <= ring_max):
                    z1 += 1
                    s1 = int( (zmin + z1)*Np*f - p1 )
                    #print([Np, zmin, z1, p1, s1])
                    if(s1+2*p1 <= ring_max): 
                            #print([Np,f,p1,s1, theta, cos(theta)])
                            s1p1 = s1 + p1
                            r1 = s1p1 + p1
                            ratio_numer = 1.0 + r1/s1
                            
                            #for z2 in range(0,zmax):
                            s2=0
                            z2=-1
                            while(s2+2*teeth_min <= ring_max):
                                    z2 += 1    
                                    for p2 in range(teeth_min, p_max):
                                        z2min=zmin_fun(p2)
                                        s2 = int( (z2min + z2)*Np*f - p2 )
                                        
                                        s2p2 = s2+p2
                                        m2 = s2p2/s1p1
                                        r2 = s2p2 + p2
                                        if m2 > m2_min and m2 < m2_max and r2 < ring_max:
                                            ratio_denom = (1- p2*r1/(p1*r2))
                                            ratio = 0 if ratio_denom == 0 else  ratio_numer/ratio_denom
                                            if ratio>ratio_min and ratio<ratio_max:
                                                results.append((Np, s1, p1, r1, s2, p2, r2, ratio, m2))
                                                #print([Np, s1, p1, r1, s2, p2, r2, ratio, m2])
                                                if len(results)>results_max:
                                                        print("Warning: max number of results exceeded")
                                                        return results
    
        return results


if __name__ == "__main__":
        results = ratio_search(teeth_min = 10, 
                ring_max = 50, 
                planets_min = 3, 
                planets_max = 6, 
                m2_min = 0.9, 
                m2_max = 1.1,
                ratio_min = 1000,
                ratio_max = 1e10)
                
        for r in results:
                print(r)
        
