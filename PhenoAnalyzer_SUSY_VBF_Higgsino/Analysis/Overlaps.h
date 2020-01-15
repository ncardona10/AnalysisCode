/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Removes particle overlaps                     
*/
#include "../ROOTFunctions.h"
#include "../DelphesFunctions.h"
#include "./Physics.h"

bool particleOverlap(vector<TLorentzVector> tlvArray)
{
  if (tlvArray.size() > 1)
  {

    bool ans = false;

    for (int i = 0; (unsigned) i < tlvArray.size() - 1 && !ans; i++)
    {
      TLorentzVector tlv1 = tlvArray[i];

      for (int j = i + 1; (unsigned) j < tlvArray.size() && !ans; j++)
      {

        TLorentzVector tlv2 = tlvArray[j];

        double deltaR = dR(tlv1, tlv2);
        ans = overlap(deltaR);
      }
    }

    return ans;
  }
  else
  {
    return false;
  }
}