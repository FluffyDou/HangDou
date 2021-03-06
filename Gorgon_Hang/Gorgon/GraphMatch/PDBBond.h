// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   A data structure to store PDB Bond information



#ifndef GRAPH_MATCH_PDB_BOND_H
#define GRAPH_MATCH_PDB_BOND_H

using namespace wustl_mm::MathTools;

namespace wustl_mm {
    namespace GraphMatch {
        class PDBBond {
        public:
            PDBBond();
            PDBBond(unsigned long long atom0Index, unsigned long long atom1Index, bool selected = false);
        public:
            unsigned long long GetAtom0Ix();
            unsigned long long GetAtom1Ix();
            bool GetSelected();
            bool GetVisible();
            void SetAtom0Ix(unsigned long long ix);
            void SetAtom1Ix(unsigned long long ix);
            void SetSelected(bool selected);
            void SetVisible(bool isVisible);
            bool tempDeleted;
            bool original;
            bool tempNew;
            bool maxOn;
            bool minOn;
            bool saddleOn;
            bool hideOn;
            float saliencies[3];
            float intensity;
            float eigenvalue;
            float eigenValue0;
            float eigenValue1;

            float logEigen0;
            float logEigen1;
            float logEigen2;

            float eigenVector0[3];
            float eigenVector1[3];
            float eigenVector2[3];

            float eigenValue0Rank;
            float eigenValue1Rank;
            float eigenValue2Rank;

            int helixIndex;
            vector<unsigned long long> helixIndices;
            bool isDisplay;
        private:
            unsigned long long atom0Ix;
            unsigned long long atom1Ix;
            bool selected;
        };

        PDBBond::PDBBond() {
            this->atom0Ix = 0;
            this->atom1Ix = 0;
            this->selected = false;
            this->tempDeleted = false;
            this->tempNew = false;
            this->original = true;
            this->maxOn = false;
            this->minOn = false;
            this->saddleOn = false;
            this->hideOn = false;
            this->helixIndex = -1;
            this->isDisplay = false;
        }

        PDBBond::PDBBond(unsigned long long atom0Index, unsigned long long atom1Index, bool selected) {
            this->atom0Ix = atom0Index;
            this->atom1Ix = atom1Index;
            this->selected = selected;
            this->tempDeleted = false;
            this->tempNew = false;
            this->original = true;
            this->maxOn = false;
            this->minOn = false;
            this->saddleOn = false;
            this->hideOn = false;
            this->helixIndex = -1;
            this->isDisplay = false;
        }


        unsigned long long PDBBond::GetAtom0Ix() {
            return atom0Ix;
        }

        unsigned long long PDBBond::GetAtom1Ix() {
            return atom1Ix;
        }

        bool PDBBond::GetSelected() {
            return selected;
        }


        void PDBBond::SetAtom0Ix(unsigned long long ix) {
            atom0Ix = ix;
        }

        void PDBBond::SetAtom1Ix(unsigned long long ix){
            atom1Ix = ix;
        }

        void PDBBond::SetSelected(bool selected) {
            this->selected = selected;
        }

    }
}
#endif
