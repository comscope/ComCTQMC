#include "Utilities.h"

namespace mpi {		
    
    Cout::Cout() : mode_(cout_mode::one) {}
    void Cout::operator=(cout_mode mode) {mode_ = mode;}
    std::ostream& Cout::operator<<(std::ostream& (*pf)(std::ostream&)) {
        if(mode_ == cout_mode::every) return pf(std::cout);
        if(mode_ == cout_mode::one) return rank() == master ? pf(std::cout) : null_;
        return null_;
    }
    
    std::ostream Cout::null_(0);
    
    Cout cout;
    
    jsx::value mpi_structure(){
        
        jsx::value jMPIStructure;
        
        std::vector<int> rankOnNode(number_of_workers());
        std::vector<int> rankOfNode(number_of_workers());
        int nNodes=0;
        
        std::vector<char> nodeNames = processor_name();
        gather(nodeNames, master);
        
        if(rank() == master) {
            
            std::map<std::string,int> rank_on_node;
            std::map<std::string,int> rank_of_node;
            
            for(int rank = 0; rank < number_of_workers(); ++rank) {
                std::string const nodeName(&nodeNames[rank*processor_name_size()], processor_name_size());
            
                if(!rank_on_node.count(nodeName)) rank_on_node[nodeName] = 0;
                else rank_on_node[nodeName]++;
                
                if(!rank_of_node.count(nodeName)) rank_of_node[nodeName] = nNodes++;
                
                rankOnNode[rank] = rank_on_node.at(nodeName);
                rankOfNode[rank] = rank_of_node.at(nodeName);
                
            }
        }
        
        all_reduce<op::sum>(rankOnNode);
        all_reduce<op::sum>(rankOfNode);
        all_reduce<op::sum>(nNodes);
        
        jMPIStructure["rank on node"] = rankOnNode[rank()];
        jMPIStructure["rank of node"] = rankOfNode[rank()];
        jMPIStructure["number of nodes"] = nNodes;
        
        return jMPIStructure;
        
    }

}
