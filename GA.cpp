/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# Implementation of the Genetic Algorithm (GA) algorithm for structure optimization
#
# Version 1.0 Octuber 2020
#
# GA algorithm has been implemented using C++;
# coupled to FHI-aims X.X (DFT code as calculator)
# It also works with newest XX version, XX
#
# Author : Jorge Refugio Fabila Fabian <jorge_fabila@ciencias.unam.mx> (IF-UNAM)
# Advisor : Dr. Oliver Paz Borbon <oliver_paz@fisica.unam.mx> (IF-UNAM)
# Advisor : Dr. Fernando Buendia Zamudio <ferbuza@fisica.unam.mx> (IF-UNAM)
#
# Note: Output folders will be generated in current directory
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/

#include"atomic.hpp"

string Simbolo_1, Simbolo_2, file_name, command, aux,geometry_file;
string initialization_file, outputfile, m_str, i_str, E_str, tag, path;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, previus;
int m, lj, N_Simbolo_1, N_Simbolo_2, count, fail_counter=0, resto, failed_max,crystal;
int n_pop, element, fit_function, init;
float step_width, Temperature, Energy, Energia, EnergiaActual, EnergiaAnterior, delta_E, k_BT, damp ;
float x_min,y_min,z_min,x_max,y_max,z_max;
Cluster clus_1, clus_2, c_aux;
Crystal cristal;
float dist, mate_mutate_ratio, swap_ratio;
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
mate_mutate_ratio=float_pipe("grep 'mate_mutate_ratio' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
swap_ratio=float_pipe("grep 'swap_ratio' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
Ncore=int_pipe("grep 'Ncore' input.bh | head -1 | cut -d \"=\" -f2 | awk '{print $1}' ");
iteraciones=int_pipe("grep 'iterations' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
swap_step=int_pipe("grep 'swap_step' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
lj=int_pipe("grep 'lennard-jones_aid' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
crystal=int_pipe("cd input ; if [ -f crystal.in ]  ; then echo 1  ;  fi ");
delta_E=float_pipe("grep 'delta_E' input.bh | cut -d \"=\" -f2 | awk '{print $1}' ");
Cluster clus[n_pop], new_cluster;
double Fit[n_pop], Probabilities[n_pop+1], rho[n_pop], a_exp, a_lin, fit_max, fit_min;
double Energies[n_pop], max_tmp, min_tmp, current, new_cluster_energy;
int seleccion[n_pop], elegido, elegido2, id_min, id_max;
float eleccion;
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
   cout<<" --> Restarting algorithm ...  "<<endl;
   string iteration_counter_i ="cd ";
   iteration_counter_i+=file_name;
   iteration_counter_i+=" ; ls Generation* | wc -l";
   i=int_pipe(iteration_counter_i,1);
   if(i==1)
   {
      string iteration_counter_m ="cd ";
      iteration_counter_m+=file_name;
      iteration_counter_m+="/Generation"+to_string(i)+" ; ls -R | grep \"relaxed.xyz\" | wc -l"; //el -R nos da libertad de escoger nombres
      m=int_pipe(iteration_counter_m,0);
      /*command.clear(); command=" cd "+file_name+"/Generation"+to_string(i)+" ; head -2 current_minimum.xyz | tail -1 | awk '{print $6 }' ";
      //RECUERDA: en cada generacion poner un current_minimum
      Energy=float_pipe(command);*/
      command.clear();
   }
   //cout<<" --> Last Generation i="<<i<<" ; last rejected m="<<m<<" ; total performed steps : "<<i+m<<endl;
   //cout<<" --> Last Generation i="<<i<<" ; "<<m<<"/"<<n_pop<<" relaxations performed : "<<i<<endl;
   cout<<" --> Last Generation i="<<i<<endl;
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
   command+=" ; fi ; mkdir "+file_name+" ; cd "+file_name+"  ; mkdir Generation1 ; cd Generation1 ; mkdir tmp_dir ; cd .. ; ";
   command+=" cp ../input/* . ;";
   system(command.c_str());
   i=1; m=0;
   contenido=0;
   init=1;
   //while(contenido!=1)
   //{
      if(initialization_file.length() > 5 ) //&& init==1 )
      {
        cout<<" --> Reading initialization file from:  "<<initialization_file<<endl;
        //Generates geometry.in and then run FHI-aims, if geometry.in.next_step is not
        //genereted then anyway here is created as a copy of the original.
        /////////////////////////////////
        command.clear(); command="cd "+file_name+"/Generation1 ; mkdir E0 ; cd ../.. ; ";
        command+="cp ";                //
        command+=initialization_file;  //
        command+=" "+file_name+"/Generation1/E0/geometry.in";
        system(command.c_str());       //
        command.clear();               //
        /////////////////////////////////
        command ="cp ";                //
        command+=initialization_file;  //
        command+=" "+file_name+"/Generation1/E0/geometry.in.next_step";
        system(command.c_str());       //
        command.clear();               //
        /////////////////////////////////
        cout<<" --> Generating a random population and adding initialization file" <<endl;
        cout<<"   --> Cluster 0 is the initializacion file from "<<initialization_file<<endl;

        count=1;
        init++;
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
               cout<<"   --> Cluster "<<element<<" created using fully random generator "<<endl;
               clus[element].srand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"     --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"   --> Cluster "<<element<<" created cleaving Au80 cluster until get the required number of atoms"<<endl;
               clus[element].rand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               clus[element].kick(step_width);
               if(lj!=0)
               {
                  cout<<"     --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"   --> Cluster "<<element<<" created using random generator based on Roy Jhonston "<<endl;
               clus[element].roy_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"     --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
         }
         else //Monometallic cases
         {
            if(randomness==1)  // fully random
            {
               cout<<"   --> Cluster "<<element<<" created using fully random generator "<<endl;
               clus[element].srand_generator(Simbolo_1,N_Simbolo_1);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"   --> Cluster "<<element<<" created cleaving Au80 cluster until get the required number of atoms"<<endl;
               clus[element].rand_generator(Simbolo_1,N_Simbolo_1);
               clus[element].kick(step_width);
               if(lj!=0)
               {
                  cout<<"   --> Optimizing geometry with L-J potential "<<endl;
                  clus[element].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"   --> Cluster "<<element<<" created using random generator based on Roy Jhonston "<<endl;
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
            geometry_file=file_name+"/Generation1/E"+to_string(element)+"/geometry.temp";
            clus[element].print_fhi(geometry_file); command.clear();
            command="cat "+file_name+"/crystal.in > "+file_name+"/Generation1/E"+to_string(element)+"/geometry.in ; ";
            command+=" sed '/atom/a initial_moment 0.5' "+geometry_file+" >> "+file_name+"/Generation1/E"+to_string(element)+"/geometry.in";
            command+=" ; rm "+geometry_file;
            system(command.c_str());
            command.clear();
         }
      }
      m=0;
}

if(i==1)
{
   cout<<" --> Starting FHI-aims calculations for first generation "<<endl;
   command="echo 'Step ----> Energy[eV]' >> "+file_name+"/Generation1/energies.txt ";
   system(command.c_str());
   command.clear();
   // running
   for(m=m;m<n_pop;m++)
   {
      cout<<"   --> Performing relaxation of element "<<m<<endl;
      command.clear();
      command="cd "+file_name+"/Generation1/E"+to_string(m)+" ; cp ../../run.sh .";
      command+=" ; cp ../../control.in .";
      system(command.c_str());
      command.clear();
      command="cd "+file_name+"/Generation1/E"+to_string(m)+" ; ./run.sh";
      system(command.c_str());
      command.clear();
      command="grep 'Have a nice day' "+file_name+"/Generation1/E"+to_string(m)+"/output.out | wc -l";
      contenido=int_pipe(command.c_str());
      command.clear();

      while(contenido!=1)
       {
          cout<<"     --> Failed SCF of element "<<m<<". Generating a new structure"<<endl;
         if(N_Simbolo_2>0)  // For bimetallic cases
         {
            if(randomness==1)  // Fully random
            {
               cout<<"     --> Cluster "<<m<<" created using fully random generator "<<endl;
               clus[m].srand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"     --> Cluster "<<m<<" created cleaving Au80 cluster until get the required number of atoms"<<endl;
               clus[m].rand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               clus[m].kick(step_width);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"     --> Cluster "<<m<<" created using random generator based on Roy Jhonston "<<endl;
               clus[m].roy_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
         }
         else //Monometallic cases
         {
            if(randomness==1)  // fully random
            {
               cout<<"     --> Cluster "<<m<<" created using fully random generator "<<endl;
               clus[m].srand_generator(Simbolo_1,N_Simbolo_1);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
            else if(randomness==0)//pseudorandomly (cuts Au80 cluster)
            {
               cout<<"     --> Cluster "<<m<<" created cleaving Au80 cluster until get the required number of atoms"<<endl;
               clus[m].rand_generator(Simbolo_1,N_Simbolo_1);
               clus[m].kick(step_width);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
            else if(randomness==2)// Roy-based generator
            {
               cout<<"     --> Cluster "<<m<<" created using random generator based on Roy Jhonston "<<endl;
               clus[m].roy_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
               if(lj!=0)
               {
                  cout<<"       --> Optimizing geometry with L-J potential "<<endl;
                  clus[m].geometry_optimization();
               }
            }
         }
         if(crystal==0)
         {
            clus[m].centroid();
            geometry_file.clear();
            geometry_file=file_name+"/Generation1/E"+to_string(m)+"/geometry.in";
            clus[m].print_fhi(geometry_file);
         }
         else{
           clus[m].centroid();
           clus[m].move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-clus[m].z_min());
           geometry_file.clear();
           geometry_file=file_name+"/Generation1/E"+to_string(m)+"/geometry.temp";
           clus[m].print_fhi(geometry_file); command.clear();
           command="cat "+file_name+"/crystal.in > "+file_name+"/Generation1/E"+to_string(m)+"/geometry.in ; ";
           command+=" sed '/atom/a initial_moment 0.5' "+geometry_file+" >> "+file_name+"/Generation1/E"+to_string(m)+"/geometry.in";
           command+=" ; rm "+geometry_file;
           system(command.c_str());
           command.clear();
         }
         command.clear();
         command="cd "+file_name+"/Generation1/E"+to_string(m)+" ; ./run.sh";
         system(command.c_str());
         command.clear();
         // Tal vez contenido ya no se use
         command="grep 'Have a nice day' "+file_name+"/Generation1/E"+to_string(m)+"/output.out | wc -l";
         contenido=int_pipe(command.c_str());
         command.clear();

       }
      cout<<"   --> Relaxation of element "<<m<<" done! "<<endl;
   }
}
   cout<<" --> Initial generation: DONE! "<<endl;
   cout<<"  "<<endl;

   cout<<"================================================================================================"<<endl;
   cout<<"GA-DFT routine starts here! "<<endl;
   cout<<"Note: "<<endl;
   cout<<"For monometallic clusters: only random xyz moves and twists will be applied "<<endl;
   cout<<"For bimetallic clusters  : atomic swap will be performed randomly with probability "<<swap_ratio<<endl;
   cout<<"================================================================================================\n"<<endl;

   while( i<iteraciones+1)
   {
      // Get energies from last iteration
      cout<<" --> Reading energies from last generation "<<endl;
      for(m=0;m<n_pop;m++)
      {
cout<<"m "<<m<<endl;
         command="grep \" | Total energy of the DFT \" "+file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/output.out | awk '{print $12}' ";
         Energies[m]=double_pipe(command.c_str());
cout<<"E "<<Energies[m]<<endl;
         m_str=to_string(m);
         E_str=string_pipe(command); //Better for Energies with all the value
         command.clear();
         command=file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/geometry.in.next_step";
         if(crystal==0) // Gas phase
         {
            clus[m].read_fhi(command);
         }
         else // Crystal
         {
            if(N_Simbolo_2>0)
            {
               clus_1=extract(command,Simbolo_1);
               clus_2=extract(command,Simbolo_2);
               clus[m]  =clus_1+clus_2;
            }
            else // for monometallics:
            {
               clus[m]=extract(command,Simbolo_1);
            }
         }
         command.clear();
         command=file_name+"/Generation"+to_string(i)+"/E"+to_string(m)+"/relaxed_coordinates.xyz";
         tag.clear();
         tag=" Iteration "+m_str+" -----> Energy = "+E_str+" eV ";
  cout<<"tag "<<tag<<endl;
         clus[m].print_xyz(command,tag);
         command.clear();
         command="echo '"+to_string(m)+" ---->' "+E_str+" >> "+file_name+"/Generation"+to_string(i)+"/energies.txt";
         system(command.c_str());
         command.clear();
      }
      cout<<" --> Sorting obtained energies "<<endl;
      // Get Maximum Energy Value
      max_tmp=Energies[0];
      for(j=1;j<n_pop;j++)
      {
         current=Energies[j];
         if( current > max_tmp )
         {
            max_tmp=current;
         }
      }
cout<<" E max "<<max_tmp<<endl;
      // Get Maximum index
      for(j=0;j<n_pop;j++)
      {
         if( Energies[j] == max_tmp )
         {
            id_max=j;
         }
      }
cout<<" E max ind "<<id_max<<endl;
      // Get Minimum Energy Value
      min_tmp=Energies[0];
      for(j=1;j<n_pop;j++)
      {
         current=Energies[j];
         if( current < min_tmp )
         {
            min_tmp=current;
         }
      }
  cout<<" E min "<<min_tmp<<endl;
  EnergiaActual=min_tmp;
      // Get Minimum index
      for(j=0;j<n_pop;j++)
      {
         if( Energies[j] == min_tmp )
         {
            id_min=j;
         }
      }
cout<<" E min ind "<<id_min<<endl;
      // Calcula rho
      cout<<" --> Normalizing energies from last generation"<<endl;
      for(j=0;j<n_pop;j++)
      {
         rho[j]=(Energies[j]-min_tmp)/(max_tmp-min_tmp);
cout<<" ----------> RHO "<<rho[j]<<endl;
      }
      // rho_max=1; rho_min=0;
      cout<<" --> Calculating Fit values with ";
      if(fit_function==0) // Exponential
      {
         cout<<" exponential function"<<endl;
      }
      else if(fit_function==1) // Linear
      {
         cout<<" linear function "<<endl;
      }
      else if(fit_function==2) // tanh
      {
         cout<<" hyperbolic tangen function "<<endl;
      }
      //Calcula fit
      for(j=0;j<n_pop;j++)
      {
         if(fit_function==0) // Exponential
         {
            Fit[j]=exp(a_exp*rho[j]);
cout<<" ----------> FIT "<<Fit[j]<<endl;
         }
         else if(fit_function==1) // Linear
         {
            Fit[j]=1-(a_lin*rho[j]);
cout<<" ----------> FIT "<<Fit[j]<<endl;
         }
         else if(fit_function==2) // tanh
         {
            Fit[j]=(0.5)*(1-tanh(2.0*rho[j]-1));
cout<<" ----------> FIT "<<Fit[j]<<endl;
         }
      }
      /*
            if((EnergiaAnterior-EnergiaActual)<delta_E)
            {
               cout<<"  --> Stopping Genetic Algorithm: delta E reached "<<endl;
               break;
            } // else:  continua con una nueva generacion
      */
      // Calcula las Probabilities
      cout<<" --> Calculating probabilities to be choosen for each element "<<endl;
      Probabilities[0]=0.0;
cout<<" ----------> PROB "<<0<<" "<<Probabilities[0]<<endl;
      for(j=1;j<n_pop+1;j++)
      {
         Probabilities[j]=Fit[j-1]+Probabilities[j-1];
cout<<" ----------> PROB "<<j<<" "<<Probabilities[j]<<endl;
      }
      cout<<" --> Chosing a random element of the population "<<endl;
      // rouleta
      eleccion=random_number(0,Probabilities[n_pop]);
cout<<" eleccion random "<<random_number(0,Probabilities[n_pop])<<"  "<<n_pop<<endl;;
cout<<eleccion <<endl;
      for(j=0;j<n_pop;j++)
      {
cout<<Probabilities[j]<<" --> "<<Probabilities[j+1]<<endl;
         if((Probabilities[j] <= eleccion) && (eleccion < Probabilities[j+1]))
         {
            elegido=j;
         }
      }
cout<<" elegido "<<elegido<<endl;
      contenido=0;
cout<<" prob max"<<Probabilities[n_pop]<<endl;
      while(contenido!=1)
      {
        cout<<"Starting new while cycle"<<endl;
         if( random_number(0,1)<mate_mutate_ratio ) //Then mate
         {
            //Code for mating
            elegido2=elegido;
            while (elegido2==elegido)
            {
               // Elige al segundo padre
               // Other rouleta
    //  cout<<" elegido2 "<<elegido2<<endl;
               eleccion=random_number(0,Probabilities[n_pop]);
  //cout<<" eleccion "<<eleccion<<endl;
               for(j=0;j<n_pop;j++)
               {
                  if(((Probabilities[j] < eleccion) && (eleccion < Probabilities[j+1]) ) )
                  {
                     elegido2=j;
  //                   cout<<"  "<<elegido2<<endl;

                  }
               }
            }
            cout<<" --> Mating choosen elements: "<<elegido<<" with "<<elegido2<<endl;
            new_cluster=Crossover(clus[elegido],clus[elegido2]);
      cout<<" performed crossover "<<endl;
            ////////////////// Bloque de código a copiar en mutate //////////////////
            path=file_name+"/Generation"+to_string(i);
            command.clear(); command="cd "+path+" ; cp ../run.sh tmp_dir/";
            command+=" ; cp ../control.in tmp_dir/";
            system(command.c_str());
            ///////////////
            if(crystal==0)
            {
               command.clear();
               command=path+"/tmp_dir/geometry.in";
               new_cluster.centroid();
               new_cluster.print_fhi(command);
               command.clear();
            }
            else
            {
              //codigo para cristal
              new_cluster.centroid();
              new_cluster.move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-new_cluster.z_min());
              geometry_file.clear(); geometry_file=file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.tmp";
              new_cluster.print_fhi(geometry_file);
              command="cat "+file_name+"/crystal.in > "+file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.in ; ";
              command+=" sed '/atom/a initial_moment 0.5' "+geometry_file+" >> "+file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.in";
              command+=" ; rm "+geometry_file;
              system(command.c_str());
              command.clear();
            }
            command.clear();
            cout<<" --> Relaxing son element "<<endl;
            command="cd "+path+"/tmp_dir ; ./run.sh";
            system(command.c_str());
            command.clear();
            command="grep 'Have a nice day' "+path+"/tmp_dir/output.out | wc -l";
            contenido=int_pipe(command.c_str());
            command.clear();
            command="grep \" | Total energy of the DFT \" "+path+"/tmp_dir/output.out | awk '{print $12}' ";
            new_cluster_energy=double_pipe(command.c_str());
            /////////////////////////////////////////////////////////////////////////
         }
         else // Assumes is mutate
         {
            if(N_Simbolo_2>0) // Si es bimetálico se utiliza el swap
            {
               if(random_number(0,1)<swap_ratio) //swap
               {
                  if(N_Simbolo_1>=N_Simbolo_2)
                  {
                     clus[elegido].swap(N_Simbolo_2);
                     cout<<" --> Swapping atoms of selected element "<<elegido<<endl;
                  }
                  else
                  {
                     clus[elegido].swap(N_Simbolo_1);
                     cout<<" --> Swapping atoms of selected element "<<elegido<<endl;
                  }
                  clus[elegido].print_xyz("tmp.xyz");
                  new_cluster.read_xyz("tmp.xyz");
                  system("rm tmp.xyz");
               }
               else //kick type
               {
                  if(random_number(0,1)<0.5) //kick
                  {
                     clus[elegido].kick(step_width);
                     clus[elegido].print_xyz("tmp.xyz");
                     new_cluster.read_xyz("tmp.xyz");
                     system("rm tmp.xyz");
                     cout<<" --> Kicking atoms of selected element "<<elegido<<endl;
                  }
                  else // twist
                  {
                     new_cluster=Crossover(clus[elegido],clus[elegido]);
                     cout<<" --> Twisting selected cluster "<<elegido<<endl;
                  }
               }
            }
            else // Significa que es monometálico
            {
               if(random_number(0,1)<0.5) //kick
               {
                  clus[elegido].kick(step_width);
                  clus[elegido].print_xyz("tmp.xyz");
                  new_cluster.read_xyz("tmp.xyz");
                  system("rm tmp.xyz");
                  cout<<" --> Kicking atoms of selected element "<<elegido<<endl;
               }
               else // twist
               {
                  new_cluster=Crossover(clus[elegido],clus[elegido]);
                  cout<<" --> Twisting selected cluster "<<elegido<<endl;
               }
            }
            ////////////////// Bloque de código a copiar en mutate //////////////////
  cout<<" empezando codigo bloque "<<endl;
            path.clear();
            path=file_name+"/Generation"+to_string(i);
            command.clear(); command="cd "+path+" ; cp ../run.sh tmp_dir/";
            command+=" ; cp ../control.in tmp_dir/";
            system(command.c_str());
            ///////////////
            if(crystal==0)
            {
               command.clear();
               command=path+"/tmp_dir/geometry.in";
               new_cluster.centroid();
               new_cluster.print_fhi(command);
               command.clear();
            }
            else
            {
              //codigo para cristal
  cout<<"codigo cristal"<<endl;
              new_cluster.centroid();
              new_cluster.move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-new_cluster.z_min());
              geometry_file.clear(); geometry_file=file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.tmp";
              new_cluster.print_fhi(geometry_file);
              command="cat "+file_name+"/crystal.in > "+file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.in ; ";
              command+=" sed '/atom/a initial_moment 0.5' "+geometry_file+" >> "+file_name+"/Generation"+to_string(i)+"/tmp_dir/geometry.in";
              command+=" ; rm "+geometry_file;
              system(command.c_str());
              command.clear();
            }
            command.clear();
            cout<<" --> Relaxing son element "<<endl;
            command="cd "+path+"/tmp_dir ; ./run.sh";
            system(command.c_str());
            command.clear();
            command="grep 'Have a nice day' "+path+"/tmp_dir/output.out | wc -l";
            contenido=int_pipe(command.c_str());
            command.clear();
            command="grep \" | Total energy of the DFT \" "+path+"/tmp_dir/output.out | awk '{print $12}' ";
            new_cluster_energy=double_pipe(command.c_str());
            /////////////////////////////////////////////////////////////////////////
      }
      if((contenido==1) && (new_cluster_energy<max_tmp))
      {
         if(i+1 <= iteraciones)
         {
            cout<<"================================================================================================"<<endl;
            cout<<" --> Starting generation "<<to_string(i+1)<<endl;
            cout<<"================================================================================================"<<endl;

         // Writing the minimum energy of the current generation
         command.clear(); command=" cp -r "+file_name+"/Generation"+to_string(i)+" "+file_name+"/Generation"+to_string(i+1);
         system(command.c_str());
         command.clear();
         command="cp "+file_name+"/Generation"+to_string(i)+"/E"+to_string(min_tmp)+"/relaxed_coordinates.xyz "+file_name;
         command+="/Generation"+to_string(i)+"/minimum_energy.xyz";
         // Obtiene la energia anterior
         EnergiaAnterior=min_tmp;
         command.clear();
         command="cp -r "+file_name+"/Generation"+to_string(i)+" "+file_name+"/Generation"+to_string(i+1);
         system(command.c_str());
         command.clear();
         command="cd "+file_name+"/Generation"+to_string(i+1)+"/E"+to_string(id_max)+" ; ";
         command+="rm * ; mv ../tmp_dir/* .";
         system(command.c_str());
         command.clear();
         command=file_name+"/Generation"+to_string(i+1)+"/E"+to_string(id_max)+"/geometry.in.next_step";
         clus[id_max].read_fhi(command);
         command.clear();
         command=file_name+"/Generation"+to_string(i+1)+"/E"+to_string(id_max)+"/relaxed_coordinates.xyz";
         tag.clear();
         tag=" Iteration "+to_string(id_max)+" -----> Energy = "+to_string(new_cluster_energy)+" eV ";
         clus[id_max].print_xyz(command,tag);
         command.clear();
         i++;
       }
      }
      else
      {
         contenido=0;
      }
   }
 }
 cout<<" --> Maximum number of generations reached ... Stopping Genetic Algorithm"<<endl;
   return 0;
}
