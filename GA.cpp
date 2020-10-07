/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# Implementation of the Genetic Algorithm (GA) algorithm for structure optimization
#
# Version X.X Octuber 2020
#
# BH algorithm has been implemented using C++;
# coupled to FHI-aims X.X (DFT code as calculator)
# It also works with newest XX version, XX
#
# Author : Jorge Refugio Fabila Fabian <jorge_fabila@ciencias.unam.mx> (IF-UNAM)
# Advisor : Dr. Oliver Paz Borbon <oliver_paz@fisica.unam.mx> (IF-UNAM)
#
# Note: Output folders will be generated in current directory
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
#include"atomic.hpp"
string Simbolo_1, Simbolo_2, file_name, command, aux,geometry_file;
string initialization_file, outputfile, m_str, i_str, E_str, tag;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, previus;
int m, lj, N_Simbolo_1, N_Simbolo_2, count, fail_counter=0, resto, failed_max,crystal;
int n_pop, element, fit_function;
float step_width, Temperature, Energy, Energia, EnergiaAnterior, k_BT, damp ;
float x_min,y_min,z_min,x_max,y_max,z_max;
Cluster clus_1, clus_2, c_aux;
Crystal cristal;
float dist;
int main(int argc, char *argv[]){
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                    Gets data from input.bh                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
Simbolo_1=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 1 ");
Simbolo_2=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 1 ");
N_Simbolo_1=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
N_Simbolo_2=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
continue_alg=int_pipe("grep 'continue' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
n_pop=int_pipe("grep 'n_pop' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
fit_function=int_pipe("grep 'fit_function' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
initialization_file=string_pipe("grep 'initialization_file' input.bh | cut -d \"=\" -f2 | awk '{print $1}'");
randomness=int_pipe("grep 'randomness' input.bh | cut -d \"=\" -f2 | awk '{print $1}'");
kick=int_pipe("grep 'kick_type' input.bh  | cut -d \"=\" -f2 | awk '{print $1}' ");
file_name=string_pipe("grep 'directory_name' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
step_width=float_pipe("grep 'step_width' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
Temperature=float_pipe("grep 'temperature_K' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
Ncore=int_pipe("grep 'Ncore' input.bh | head -1 | cut -d \"=\" -f2 | awk '{print $1}' ");
iteraciones=int_pipe("grep 'iterations' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
swap_step=int_pipe("grep 'swap_step' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
lj=int_pipe("grep 'lennard-jones_aid' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
crystal=int_pipe("cd input ; if [ -f crystal.in ]  ; then echo 1  ;  fi ");
Cluster clus[n_pop];
double Fit[n_pop], rho[n_pop], a_exp, a_lin, fit_max, fit_min;
double Energies[n_pop], max_tmp, min_tmp, current;
int seleccion[n_pop];
// Meta-parámetros /////
failed_max=3;         //
damp=0.7;             //
dist=1.0;             //
a_exp=-3.0;           //
a_lin=0.7;            //
////////////////////////
srand(time(NULL)); // init  Randomness
//Automatically detects if exists a crystal file
if(crystal==1)  //Esto sustituye tener que poner [x_min,x_max]; [y_min,y_max]... en el input
{
  cout<<" --> Reading crystal from file "<<endl;
  cristal.read_fhi("input/crystal.in");
  x_min=cristal.x_min();
  x_max=cristal.x_max();
  y_min=cristal.y_min();
  y_max=cristal.y_max();
  z_min=cristal.z_min();
  z_max=cristal.z_max();
}
else
{
  cout<<" --> crystal.in file not found ... performing gas phase search "<<endl;
}
int i = 1;

if(continue_alg==1)
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                      RESTART ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //En vez de reinicio podriamos hacer que directamente lo genera desde aca
  cout<<" --> Restarting algorithm ...  "<<endl;
      string iteration_counter_i ="cd ";
             iteration_counter_i+=file_name;
             iteration_counter_i+=" ; ls Generation* | wc -l";
  i=int_pipe(iteration_counter_i,1);
      string iteration_counter_m ="cd ";
             iteration_counter_m+=file_name;
             iteration_counter_m+="/Generation"+to_string(i)+" ; ls -R | grep \"relaxed.xyz\" | wc -l"; //el -R nos da libertad de escoger nombres
  m=int_pipe(iteration_counter_m,0);
     command.clear(); command=" cd "+file_name+"/Generation"+to_string(i)+" ; head -2 current_minimum.xyz | tail -1 | awk '{print $6 }' ";
     //RECUERDA: en cada generacion poner un current_minimum
     Energy=float_pipe(command);
     command.clear();

  //cout<<" --> Last Generation i="<<i<<" ; last rejected m="<<m<<" ; total performed steps : "<<i+m<<endl;
  cout<<" --> Last Generation i="<<i<<" ; "<<m<<"/"<<n_pop<<" relaxations performed : "<<i<<endl;
  cout<<" --> Restarting from generation"<<i<<" and relaxation "<<m<<" "<<endl;
//  i++;
//  cout<<" --> Starting step "<<i<<endl;
  m++;
}
else
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                        BEGIN ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
   cout<<" --> Starting a new search "<<endl;
   // Creates work directory
   command ="if [ -d "+file_name+" ] ; then mv "+file_name+" other_"+file_name;
   command+=" ; fi ; mkdir "+file_name+" ; cd "+file_name+"  ; mkdir Generation1 ;";
   command+=" cp ../input/* . ;";
   system(command.c_str());
   i=1; m=0;
   contenido=0;
//   while(contenido!=1)
//   {
      if(initialization_file.length() > 5)
      {
        cout<<" --> Reading initialization file from:  "<<initialization_file<<endl;
        //Generates geometry.in and then run FHI-aims, if geometry.in.next_step is not
        //genereted then anyway here is created as a copy of the original.
        /////////////////////////////////
        command.clear(); command="cd "+file_name+"/Generation1 ; mkdir E"+to_string(m);
        command+="cp ";                //
        command+=initialization_file;  //
        command+=" "+file_name+"/Generation1/E1/geometry.in";
        system(command.c_str());       //
        command.clear();               //
        /////////////////////////////////
        command ="cp ";                //
        command+=initialization_file;  //
        command+=" "+file_name+"/Generation1/E1/geometry.in.next_step";
        system(command.c_str());       //
        command.clear();               //
        /////////////////////////////////
        cout<<" --> Generating a random population and adding initialization file from "<<initialization_file<<endl;
        count=1;
      }
      else
      {
//Genera Generation1/E1/geometry.in aleatorio
         cout<<" --> Generating a random population "<<endl;
         count=0;
      }
      for(element=count;element<n_pop;element++)
      {
         if(N_Simbolo_2>0)  // For bimetallic cases
         {
            if(randomness==1)  // Fully random
            {
               cout<<"   --> Using fully random generator "<<endl;
               clus[element].srand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"   --> Cleaving Au80 cluster until get a "<<Simbolo_1<<N_Simbolo_1<<Simbolo_2<<N_Simbolo_2<<" new cluster "<<endl;
               clus[element].rand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"   --> Using random generator based on Roy Jhonston "<<endl;
               clus[element].roy_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
         }
         else //Monometallic cases
         {
            if(randomness==1)  // fully random
            {
               cout<<"   --> Using fully random generator "<<endl;
               clus[element].srand_generator(Simbolo_1,N_Simbolo_1);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"   --> Cleaving Au80 cluster until get a "<<Simbolo_1<<N_Simbolo_1<<" new cluster "<<endl;
               clus[element].rand_generator(Simbolo_1,N_Simbolo_1);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"   --> Using random generator based on Roy Jhonston "<<endl;
               clus[element].roy_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
         }
         command.clear(); command="cd "+file_name+"/Generation1 ; mkdir E"+to_string(element);
         system(command.c_str());
         command.clear();
         if(crystal==0)
         {
            clus[element].centroid();
            geometry_file.clear();
            geometry_file=file_name+"/Generation1/E"+to_string(element)+"/geometry.in";
            clus[element].print_fhi(geometry_file);
         }
         else{
            clus[element].centroid();
            clus[element].move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-clus[element].z_min());
            geometry_file.clear();
            geometry_file=file_name+"/Generation1/E"+to_string(element)+"/geometry.in";
            clus[element].print_fhi(geometry_file);
            command="cat "+file_name+"/crystal.in > "+geometry_file+" ; sed '/atom/a initial_moment 0.5' "+file_name;
            command+="/geometry.tmp >> "+geometry_file+" ; rm "+file_name+"/geometry.tmp" ;
            system(command.c_str());
            command.clear();
         }
      }
      m=0;
}
      cout<<" --> Starting FHI-aims calculations "<<endl;

      command="echo 'Step ----> Energy[eV]' >> "+file_name+"/Generation1/energies.txt ";
      system(command.c_str());
      command.clear();
cout<<"Entering loop"<<endl;
// running
for(m=m;m<n_pop;m++)
{
  cout<<m<<endl;
      command.clear();
      command="cd "+file_name+"/Generation1/E"+to_string(m)+" ; cp ../../run.sh .";
      command+=" ; cp ../../control.in .";
      system(command.c_str());
      command.clear();
      command="cd "+file_name+"/Generation1/E"+to_string(m)+" ; ./run.sh";
      system(command.c_str());
      command.clear();
      // Tal vez contenido ya no se use
      command="grep 'Have a nice day' "+file_name+"/Generation1/E"+to_string(m)+"/output.out | wc -l";
      contenido=int_pipe(command.c_str());
      command.clear();
}

   cout<<" --> Initial generation: DONE! "<<endl;
   cout<<" --> AG-DFT routine starts here "<<endl;
   cout<<"  "<<endl;
   cout<<" --> Starting generation 2 "<<endl;

   cout<<"================================================================================================"<<endl;
   cout<<"GA-DFT routine starts here! "<<endl;
   /*cout<<"Note: "<<endl;
   cout<<"For monometallic clusters: only random xyz moves will be applied "<<endl;
   cout<<"For bimetallic clusters  : 1 atomic swap will be performed after "<<swap_step<<" moves "<<endl;
   cout<<"================================================================================================"<<endl;*/

   while( i<iteraciones)
   {
      // Get energies from last iteration
      for(m=0;m<n_pop;m++)
      {
         command="grep \" | Total energy of the DFT \" "+file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/output.out | awk '{print $12}' ";
         Energies[m]=double_pipe(command.c_str());
         m_str=to_string(m);
         E_str=string_pipe(command); //Better for Energies with all the value
         command.clear();
         command=file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/geometry.in.next_step";
         clus[element].read_fhi(command);
         command.clear();
         command=file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/relaxed_coordinates.xyz";
         tag.clear();
         tag=" Iteration "+m_str+" -----> Energy = "+E_str+" eV ";
         clus[element].print_xyz(command,tag);
         command.clear();
         command="echo '"+to_string(m)+" ---->' "+E_str+" >> "+file_name+"/Generation"+to_string(i)+"/energies.txt";
         system(command.c_str());
         command.clear();
      }
      // Sort energies
         // Get Maximum
      max_tmp=Energies[0];
      for(j=1;j<n_pop;j++)
      {
         current=Energies[j];
         if( current > max_tmp )
         {
            max_tmp=current;
         }
      }
         // Get Minimum
      min_tmp=Energies[0];
      for(j=1;j<n_pop;j++)
      {
         current=Energies[j];
         if( current > min_tmp )
         {
            min_tmp=current;
         }
      }
      // Calcula rho
      for(j=1;j<n_pop;j++)
      {
         rho[j]=(Energies[j]-min_tmp)/(max_tmp-min_tmp);
      }
      // rho_max=1; rho_min=0;

      //Calcula fit
      for(j=1;j<n_pop;j++)
      {
         if(fit_function==0) // Exponential
         {
            Fit[j]=exp(a_exp*rho[j]);
         }
         else if(fit_function==1) // Linear
         {
            Fit[j]=1-(a_lin*rho[j]);
         }
         else if(fit_function==2) // tanh
         {
            Fit[j]=(0.5)*(1-tanh(2.0*rho[j]-1));
         }
      }
/*
      if(criteria==0)
      {
         break;
      }
*/
      //Establece los parámetros, rangos y porcentajes
      if(fit_function==0) // Exponential
      {
         fit_max=exp(a_exp*1.0);
         fit_min=exp(a_exp*0.0);
      }
      else if(fit_function==1) // Linear
      {
         fit_max=1-(a_lin*1.0);
         fit_min=1-(a_lin*0.0);
      }
      else if(fit_function==2) // tanh
      {
         fit_max=(0.5)*(1-tanh(2.0-1));
         fit_min=(0.5)*(1-tanh(-1.0));
      }

      //Seleccion de los especímenes
      for(j=1;j<n_pop;j++)
      {
         if(Fit[j]>(fit_min+(fit_max-fit_min)*0.8))  // 20% sobrevive tal cual
         {
            seleccion[j]=0; //"sobrevive"
         }
         else if((Fit[j]>(fit_min+(fit_max-fit_min)*0.5)) && (Fit[j]<(fit_min+(fit_max-fit_min)*0.8))) // %30 pasa sus hijos
         {
            seleccion[j]=1; //"crossover"
         }
         else if(Fit[j]<=(fit_min+(fit_max-fit_min)*0.5)) // %50 mutan: 20% kick; 10% swap; 20% twist
         {
            if(N_Simbolo_2>0) // Si es bimetálico se utiliza el swap
            {
               if(random_number(0,1)<0.75)
               {
                  seleccion[j]=2; //"kick";
               }
               else
               {
                  seleccion[j]=3; //"swap";
               }
            }
            else
            {
               seleccion[j]=2; //"kick";
            }
         }
      }
      // Performs opereations
      for(j=1;j<n_pop;j++)
      {
         if(seleccion[j]==0)  // 20% sobrevive tal cual
         {

         }
         else if(seleccion[j]==1) // %30 pasa sus hijos
         {

         }
         else if(seleccion[j]==2) // %50 mutan: 20% kick; 10% swap; 20% twist??
         {
         }
         else if(seleccion[j]==3) // swap
         {

         }
      }
      // Relaja las configurasaoes
      for(j=1;j<n_pop;j++)
      {
      }
      i++;
   }
   return 0;
}
