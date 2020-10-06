#include"atomic.hpp"

Cluster clus1, clus2, clus3;
int main()
{
clus1.srand_generator("Au",500);
clus2.srand_generator("Ag",500);
clus3=Crossover(clus1,clus2);
clus3.print_xyz("clus3.xyz");
//clus3.show("VESTA");

return 0;
}
