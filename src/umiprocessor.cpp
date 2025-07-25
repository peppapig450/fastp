#include "umiprocessor.h"

UmiProcessor::UmiProcessor(Options* opt){
    mOptions = opt;
}


UmiProcessor::~UmiProcessor(){
}

void UmiProcessor::process(Read* r1, Read* r2) {
    if(!mOptions->umi.enabled)
        return;

    const auto& r1Len       = static_cast<int>(r1->length());
    const auto& r2Len       = static_cast<int>(r2->length());
    const auto& umiLocation = mOptions->umi.location;
    const auto& umiLen      = mOptions->umi.length;
    const auto& umiSkip     = mOptions->umi.skip;

    string umi;
    if(umiLocation == UMI_LOC_INDEX1)
        umi = r1->firstIndex();
    else if(umiLocation == UMI_LOC_INDEX2 && r2)
        umi = r2->lastIndex();
    else if(umiLocation == UMI_LOC_READ1){
        umi = r1->seq().substr(0, std::min(r1Len, umiLen));
        r1->trimFront(umi.length() + umiSkip);
    }
    else if(umiLocation == UMI_LOC_READ2 && r2){
        umi = r2->seq().substr(0, std::min(r2Len, umiLen));
        r2->trimFront(umi.length() + umiSkip);
    }
    else if(umiLocation == UMI_LOC_PER_INDEX){
        string umiMerged = r1->firstIndex();
        if(r2) {
            umiMerged = umiMerged + "_" + r2->lastIndex();
        }

        addUmiToName(r1, umiMerged);
        if(r2) {
            addUmiToName(r2, umiMerged);
        }
    }
    else if(umiLocation == UMI_LOC_PER_READ){
        string umiMerged = r1->seq().substr(0, std::min(r1Len, umiLen));
        r1->trimFront(umiMerged.length() + umiSkip);
        if(r2){
            string umi2 = r2->seq().substr(0, std::min(r2Len, umiLen));
            umiMerged = umiMerged + "_" + umi2;
            r2->trimFront(umi2.length() + umiSkip);
        }

        addUmiToName(r1, umiMerged);
        if(r2){
            addUmiToName(r2, umiMerged);
        }
    }

    if(umiLocation != UMI_LOC_PER_INDEX && umiLocation != UMI_LOC_PER_READ) {
        if(r1 && !umi.empty()) 
            addUmiToName(r1, umi);
        if(r2 && !umi.empty())
            addUmiToName(r2, umi);
    }
}

void UmiProcessor::addUmiToName(Read* r, string umi){
    string tag;
    string delimiter = mOptions->umi.delimiter;
    if(mOptions->umi.prefix.empty())
        tag = delimiter + umi;
    else
        tag = delimiter + mOptions->umi.prefix + "_" + umi;
    int spacePos = -1;
    for(int i=0; i<r->mName->length(); i++) {
        if(r->mName->at(i) == ' ') {
            spacePos = i;
            break;
        }
    }
    if(spacePos == -1) {
        r->mName->append(tag);
    } else {
        r->mName->insert(spacePos, tag);
    }

}


bool UmiProcessor::test() {
    return true;
}
