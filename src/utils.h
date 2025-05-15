#ifndef utils_h
#define utils_h

#include "dataset.h"

class Utils {

public:

    Utils() = default;
    ~Utils() = default;

    vector<pieceOfWork> divideWork(double numItems, int& numProcessors) {
        // divide work between processors
        vector<pieceOfWork> work;

        if (numItems < numProcessors) { numProcessors = numItems; }
        size_t startIndex = 0;

        for (size_t remainingProcessors = numProcessors; remainingProcessors > 0;
        remainingProcessors--) {

            //case for last processor
            size_t numToProcess = numItems;
            if (remainingProcessors != 1) {
                numToProcess = ceil(numItems / remainingProcessors);
            }
            work.push_back(pieceOfWork(startIndex, (startIndex+numToProcess)));
            startIndex += numToProcess;
            numItems -= numToProcess;
        }
        return work;
    }

    float round2Places(float var) {
        // 37.66666 * 100 =3766.66
        // 3766.66 + .5 =3767.16    for rounding off value
        // then type cast to int so value is 3767
        // then divided by 100 so the value converted into 37.67
        float value = (int)(var * 100 + .5);
        return (float)value / 100;
    }


private:


};

#endif /* utils_h */

