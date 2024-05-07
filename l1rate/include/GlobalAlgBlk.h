#ifndef DataFormats_L1TGlobal_GlobalAlgBlk_h
#define DataFormats_L1TGlobal_GlobalAlgBlk_h

/**
* \class GlobalAlgBlk
*
*
* Description: L1 micro Global Trigger - Block holding Algorithm Information
*
* Implementation:
* <TODO: enter implementation details>
*
* \author: Brian Winer - Ohio State
*
*
*/

// system include files
#include <vector>
#include <iostream>
#include <iomanip>

// user include files

// forward declarations


class GlobalAlgBlk;

// class interface

class GlobalAlgBlk
{
public:
    /// constructors
    GlobalAlgBlk(); // empty constructor, all members set to zero;

    GlobalAlgBlk(int orbitNr, int bxNr, int bxInEvent);

    /// destructor
    virtual ~GlobalAlgBlk();

public:
    static constexpr unsigned int maxPhysicsTriggers = 512;

    /// set simple members
    void setL1MenuUUID(int uuid)      { m_orbitNr       = uuid; }
    void setL1FirmwareUUID(int fuuid) { m_bxNr   = fuuid; }
    void setbxInEventNr(int bxNr)   { m_bxInEvent      = bxNr; }
    void setFinalORVeto(bool fOR)   { m_finalORVeto    = fOR; }
    void setFinalORPreVeto(bool fOR){ m_finalORPreVeto = fOR; }
    void setFinalOR(bool fOR)       { m_finalOR        = fOR; }
    void setPreScColumn(int psC)    { m_preScColumn    = psC; }

    /// get simple members
    inline int getL1MenuUUID() const      { return m_orbitNr; }
    inline int getL1FirmwareUUID() const  { return m_bxNr; }
    inline int getbxInEventNr() const     { return m_bxInEvent; }
    inline bool getFinalOR() const        { return m_finalOR; }
    inline bool getFinalORPreVeto() const { return m_finalORPreVeto; };
    inline bool getFinalORVeto() const    { return m_finalORVeto; }
    inline int getPreScColumn() const     { return m_preScColumn; }

    /// Copy vectors words
    void copyInitialToInterm() { m_algoDecisionPreScaled   = m_algoDecisionInitial; }
    void copyIntermToFinal() { m_algoDecisionFinal   = m_algoDecisionPreScaled; }

    /// Set decision bits
    void setAlgoDecisionInitial(unsigned int bit, bool val);
    void setAlgoDecisionInterm(unsigned int bit, bool val);
    void setAlgoDecisionFinal(unsigned int bit, bool val);

    /// Get decision bits
    std::vector<bool> const & getAlgoDecisionInitial() const {
        return m_algoDecisionInitial;
    }
    std::vector<bool> const & getAlgoDecisionInterm() const {
        return m_algoDecisionPreScaled;
    }
    std::vector<bool> const & getAlgoDecisionFinal() const { 
        return m_algoDecisionFinal;
    }
    bool getAlgoDecisionInitial(unsigned int bit) const;
    bool getAlgoDecisionInterm(unsigned int bit) const;
    bool getAlgoDecisionFinal(unsigned int bit) const;

    /// reset the content of a GlobalAlgBlk
    void reset();

    // compare the content of this GlobalAlgBlk with another one
    virtual bool operator==(const GlobalAlgBlk& rhs) const;
    virtual inline bool operator!=(const GlobalAlgBlk& rhs) const { return !(operator==(rhs)); };

    /// pretty print the content of a GlobalAlgBlk
    void print(std::ostream& myCout) const;

private:

    // where noted member data has been re-interpreted, to keep persistant data the same, as these features were added late in release cycle.

    /// orbit number -> L1MenuUUID
    int m_orbitNr;

    /// bunch cross number of the actual bx -> L1FirmwareUUID
    int m_bxNr;

    /// bunch cross in the GT event record (E,F,0,1,2)
    int m_bxInEvent;

    // finalOR 
    bool m_finalOR;
    bool m_finalORPreVeto;
    bool m_finalORVeto; 
    
    //Prescale Column
    int m_preScColumn;

    std::vector<bool> m_algoDecisionInitial;
    std::vector<bool> m_algoDecisionPreScaled;   // -> Interm
    std::vector<bool> m_algoDecisionFinal;
};

#endif /*DataFormats_L1TGlobal_GlobalAlgBlk_h*/
