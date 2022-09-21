#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

std::vector<int> coreEU0;
std::vector<int> coreEU1;
std::vector<int> coreEU2;
std::vector<int> coreEU3;
std::vector<int> corePU4;
std::vector<int> corePU5;
std::vector<int> corePU6;
std::vector<int> corePU7;

std::vector<int> clusterE;
std::vector<int> clusterP;

std::vector<int> header;

long int avgE = 0;
long int avgP = 0;
long int countE = 0;
long int countP = 0;

std::vector<int> createHeader(unsigned long size){
    std::vector<int> ret;
    for (int i = 1; i<=size; ++i) ret.push_back(i);
    return ret;
}

void writeDataToFile(std::vector<int> vector, std::string name){
    //open file for writing
    std::ofstream fw("/Users/javi/Desktop/ResultadosCores/"+name, std::ofstream::out);
    if (fw.is_open()){
      //store array contents to text file
      for (int i = 0; i < vector.size(); i++) fw << vector[i]<<"\t";
      fw.close();
    }
    else std::cout << "Problem opening file";
}

int main(){
    
    //Data entry
    std::ifstream in("datos.txt");
    auto cinbuf = std::cin.rdbuf(in.rdbuf());
    //Beginning
    std::string linea;
    while (getline(std::cin, linea)) {
        std::stringstream ss(linea);
        std::string temp;
        ss>>temp;
        //Analisis of the cpu active residencies
        
        if (temp=="cpu"){
            ss>>temp;
            if (temp=="0"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp; //residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    coreEU0.push_back(std::stoi(temp));
                }
            }
            if (temp=="1"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    coreEU1.push_back(std::stoi(temp));
                }
            }
            if (temp=="2"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    coreEU2.push_back(std::stoi(temp));
                }
            }
            if (temp=="3"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    coreEU3.push_back(std::stoi(temp));
                }
            }
            if (temp=="4"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    corePU4.push_back(std::stoi(temp));
                }
            }
            if (temp=="5"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    corePU5.push_back(std::stoi(temp));
                }
            }
            if (temp=="6"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp;//residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    corePU6.push_back(std::stoi(temp));
                }
            }
            if (temp=="7"){
                ss>>temp;
                if (temp=="active"){
                    ss>>temp; //residency:
                    ss>>temp; //3.06%
                    temp.pop_back();
                    corePU7.push_back(std::stoi(temp));
                }
            }
        }
         
        //Power consumption
        if (temp=="E-Cluster"){
            ss>>temp;
            if (temp=="Power:"){
                int consumo;
                ss>>consumo;
                avgE += consumo;
                clusterE.push_back(consumo);
                ++countE;
            }
        }
        
        if (temp=="P-Cluster"){
            ss>>temp;
            if (temp=="Power:"){
                int consumo;
                ss>>consumo;
                avgP += consumo;
                clusterP.push_back(consumo);
                ++countP;
            }
        }

    }
    
    //Done
    std::cout<<"Times counted: "<<coreEU0.size()<<"\n";
    std::cout<<"Average Power E: "<<avgE/countE<<std::endl;
    std::cout<<"Average Power P: "<<avgP/countP<<std::endl;
    std::cout<<"Total consumption: "<<avgE/countE+avgP/countP<<std::endl;
    header = createHeader(coreEU0.size());
    
    
    writeDataToFile(coreEU0, "coreEU0.txt");
    writeDataToFile(coreEU1, "coreEU1.txt");
    writeDataToFile(coreEU2, "coreEU2.txt");
    writeDataToFile(coreEU3, "coreEU3.txt");
    writeDataToFile(corePU4, "corePU4.txt");
    writeDataToFile(corePU5, "corePU5.txt");
    writeDataToFile(corePU6, "corePU6.txt");
    writeDataToFile(corePU7, "corePU7.txt");
    writeDataToFile(clusterE, "clusterE.txt");
    writeDataToFile(clusterP, "clusterP.txt");
    writeDataToFile(header, "header.txt");
    
    
    std::cin.rdbuf(cinbuf);
}
