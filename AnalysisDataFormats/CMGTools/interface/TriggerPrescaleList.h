#ifndef _AnalysisDataFormats_CMGTools_TriggerPrescaleList_H_
#define _AnalysisDataFormats_CMGTools_TriggerPrescaleList_H_

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/Ref.h"

namespace cmg {

class TriggerPrescaleList {
    public:
        TriggerPrescaleList() : triggerNames_(0) {}
        TriggerPrescaleList(const edm::Handle<edm::TriggerResults> & handle) ;
        ~TriggerPrescaleList() {}

        // get prescale by index.
        int getPrescaleForIndex(int index) const ;
        // get prescale by name or name prefix (if setTriggerNames was called)
        int getPrescaleForName(const std::string & name, bool prefixOnly=false) const ; 

        // return the TriggerResults associated with this
        const edm::TriggerResults & triggerResults() const {
            return *edm::getProduct<edm::TriggerResults>(triggerResults_);
        }

        // use this method first if you want to be able to access the prescales by name
        // you can get the TriggerNames from the TriggerResults and the Event (edm or fwlite)
        void setTriggerNames(const edm::TriggerNames &names) const { triggerNames_ = &names; }

        // set that the trigger of given index has a given prescale
        void addPrescaledTrigger(int index, int prescale) ;

    protected:
        std::vector<int> prescaleValues_;
        edm::RefCore triggerResults_;
        mutable const edm::TriggerNames *triggerNames_;
};


}

#endif
