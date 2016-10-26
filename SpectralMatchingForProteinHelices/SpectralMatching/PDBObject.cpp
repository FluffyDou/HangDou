/*********************************************************************/
// PDBObject.cpp
// This is a class holding all the data from standard PDB file.

// A PDB object has a parser to read in standard PDB file (.pdb).
// A PDB object holds a list of Chain Object, ProteinChain.

// A ProteinChain has information as that in Mathematica, like atom 
// coordinates, sequence, secondary structure sequence and so on

// Hang 10/06/2014
/*********************************************************************/

#include "PDBObject.h"
#include "Atom.h"

// Create an empty object
PDBObject::PDBObject(void)
{
    Clear();
}

// Read in all the chains from the input .pdb file.
// The input file should be standard PDB file with valid header
// The input file is allowed to be truncated, with less chains than the header claims
PDBObject::PDBObject(char* fileName)
{
    ReadInPDBFile(fileName);
}

// Read in specific chain information from the input .pdb file.
PDBObject::PDBObject(char* fileName, int chainID)
{

}

PDBObject::~PDBObject(void)
{

}

// Read in the PDB file. This operation will evoke clearance of existing chains
void PDBObject::ReadInPDBFile(const char* fileName)
{
    // Open the file
    std::ifstream file(fileName);
    
    if( !file.is_open() )
        std::cout<<"Unable to open "<<fileName<<"...\n";
    else
    {
        /****** Clear the existing chains ******/
        Clear();

        /******** Read in new stuff ********/

        // Read in line by line
        std::string line;
        
        //char  count       = 0;   // The number of chains
        //char  oldChainID  = 'Z'; // Use to control when to update residue index offset for each chain
        int   oldIndex    = 0;     // Use to control when do we start reading next residue
        bool  finishCount = false; // Finish counting chains or not


        // Set chain ID offset
        //m_chainIdOffset = 'A';   // By default, the chains start from 'A' in PDB file
        while( std::getline(file, line) )
        {
            // Parse the current line by space
            std::istringstream buf(line);
            std::istream_iterator<std::string> beg(buf), end;
            std::vector<std::string> tokens(beg, end);


            /******** Count the number of chains from the header. ********/
            // We can hold less chain than this count, but not more.
            if( "COMPND" == tokens[0] && "CHAIN:" == tokens[2] )
            {             
                //count += tokens.size() - 3;
                // The chain offset is the first chain list in "COMPND CHAIN:" row
                //m_chainIdOffset = tokens[3].c_str()[0];

                // Create empty chains
                for(auto it = tokens.begin() + 3; it != tokens.end(); ++it)
                {
                    m_proteinChains[it->c_str()[0]] = ProteinChain();
                }
            }
            /************** Create empty ProteinChain objects ************/
            //else if( "SOURCE" == tokens[0] && false == finishCount )
            //{
            //    // Since we always read in standard PDB file, when we start on "SOURCE" line,
            //    // it means we have finished counting possible chains. So we initialize chain objects                
            //    // Some ProteinChain objects may stay empty if specific chains are required
            //    // or the file only has the header but not full chain information.
            //    // So far, we don't actually read any information from SOURCE.

            //    std::vector<ProteinChain> proteinChains(count);
            //    m_chains = proteinChains;

            //    // Assign Chain IDs --- start may not start from A. They increase continuously.
            //    for(char i = 0; i < count; ++i)
            //        m_chains[i].SetChainID(m_chainIdOffset + i);

            //    finishCount = true;
            //}
             /******** Actually read in data ********/
            else if( "HELIX" == tokens[0] )
            {
                // Get the chain id, start and end residue index for this helix
                //int chainID       = readString(line,19,19).c_str()[0] - m_chainIdOffset;
                char chainID     = readString(line,19,19).c_str()[0];
                int startResIndex = readInt(line,21,24);
                int endResIndex   = readInt(line,33,36);
                // The helix really only holds the start and end residue index.
                std::vector<int> residueIndex(2);
                residueIndex[0] = startResIndex;
                residueIndex[1] = endResIndex;
                //m_chains[chainID].AddHelixToSSE( AlphaHelix(residueIndex) );

                // The protection here is actually not necessary, since a valid file header is required
                //if(m_proteinChains.find(chainID2) != m_proteinChains.end())
                m_proteinChains[chainID].AddHelixToSSE( AlphaHelix(residueIndex) );
            }
            else if( "ATOM" == tokens[0] )
            {
                // Get the chain id
                //int chainID      = readString(line,21,21).c_str()[0] - m_chainIdOffset;
                char chainID    = readString(line,21,21).c_str()[0];
                int residueIndex = readInt(line,22,25);

                // Read in sequence --- ATTENTION! index in C++ starts from 0 and index in PDB file starts from 1
                // As Mathematica, we only read in those sequences that have atom information
                if( residueIndex != oldIndex )
                {
                    char residueChar = PDBObject::GetAminoAbbrev(readString(line,17,19).c_str());                
                    //m_chains[chainID].AddSeqChar( residueChar );
                    //m_chains[chainID].AddResidue( residueIndex, Residue(residueChar) );

                    // The protection here is actually not necessary, since a valid file header is required
                    //if(m_proteinChains.find(chainID2) != m_proteinChains.end())
                    m_proteinChains[chainID].AddSeqChar( residueChar );
                    m_proteinChains[chainID].AddResidue( residueIndex, Residue(residueChar) );

                }
            
                // Add a flag --- for those empty chains not mentioned in "HELIX" rows, do not bother inserting atoms. ToDo.
                // Read the atom into corresponding residue
                Atom atom( readFloat(line,30,37),  readFloat(line,38,45), readFloat(line,46,53) );
                //m_chains[chainID].AddAtomToResidue(residueIndex, atom);

                // The protection here is actually not necessary, since a valid file header is required
                //if(m_proteinChains.find(chainID2) != m_proteinChains.end())
                m_proteinChains[chainID].AddAtomToResidue(residueIndex, atom);

                //oldChainID = chainID;
                oldIndex   = residueIndex;
            }

        }
        // Close the file
        file.close();
        
    }

}

// Retrieve information from the read in line from PDB file
std::string PDBObject::trimString(std::string const& str)
{
    std::size_t first = str.find_first_not_of(' ');    
    
    if( first == std::string::npos )
    {
        std::cout<<"Try to trim a space only string...";
        // Return an empty string
        return "";
    }

    std::size_t last  = str.find_last_not_of(' ');

    return str.substr(first, last-first+1);
}
std::string PDBObject::readString(std::string line, int start, int end)
{
    return trimString( line.substr(start, end-start+1) );
}
int PDBObject::readInt(std::string line, int start, int end)
{
    return std::stoi( line.substr(start, end-start+1) );
}
float PDBObject::readFloat(std::string line, int start, int end)
{
    return std::stof(line.substr(start, end-start+1));
}

// Get the abbreviation of an amino acid
char PDBObject::GetAminoAbbrev(const char *aminoAcid) {
    char result;
    if(strcmp(aminoAcid, "ALA") == 0) {
        result = 'A';
    } else if(strcmp(aminoAcid, "ARG") == 0) {
        result = 'R';
    } else if(strcmp(aminoAcid, "ASN") == 0) {
        result = 'N';
    } else if(strcmp(aminoAcid, "ASP") == 0) {
        result = 'D';
    } else if(strcmp(aminoAcid, "CYS") == 0) {
        result = 'C';
    } else if(strcmp(aminoAcid, "GLN") == 0) {
        result = 'Q';
    } else if(strcmp(aminoAcid, "GLU") == 0) {
        result = 'E';
    } else if(strcmp(aminoAcid, "GLY") == 0) {
        result = 'G';
    } else if(strcmp(aminoAcid, "HIS") == 0) {
        result = 'H';
    } else if(strcmp(aminoAcid, "ILE") == 0) {
        result = 'I';
    } else if(strcmp(aminoAcid, "LEU") == 0) {
        result = 'L';
    } else if(strcmp(aminoAcid, "LYS") == 0) {
        result = 'K';
    } else if(strcmp(aminoAcid, "MET") == 0) {
        result = 'M';
    } else if(strcmp(aminoAcid, "PHE") == 0) {
        result = 'F';
    } else if(strcmp(aminoAcid, "PRO") == 0) {
        result = 'P';
    } else if(strcmp(aminoAcid, "SER") == 0) {
        result = 'S';
    } else if(strcmp(aminoAcid, "THR") == 0) {
        result = 'T';
    } else if(strcmp(aminoAcid, "TRP") == 0) {
        result = 'W';
    } else if(strcmp(aminoAcid, "TYR") == 0) {
        result = 'Y';
    } else if(strcmp(aminoAcid, "VAL") == 0) {
        result = 'V';
    } else {
        printf("Invalid amnio acid name.../n");
        result = 'Z';
    }
    return result;
}
