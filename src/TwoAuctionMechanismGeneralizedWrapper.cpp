#include "TwoAuctionMechanismGeneralized.h"

#ifndef TWO_AUCTION_GENERILIZED_WRAPPER_H
#define TWO_AUCTION_GENERILIZED_WRAPPER_H

using namespace two_auction_mechanism;

extern "C"
{

	TwoAuctionMechanismGeneralized * 
	two_auction_mechanism_generalized_new()
	{
		return new TwoAuctionMechanismGeneralized();
	}
	
	int 
	two_auction_mechanism_generalized_zero_in(TwoAuctionMechanismGeneralized *auction_mechanism, 
		int nh, int nl, double bh, double bl, double kh, double kl, 
		double rph, double rpl, double Q,  double *x, double *y)
	{
		return auction_mechanism->zeroin(nh, nl, bh, bl, kh, kl, rph, rpl, Q, x, y);
	}
	
	void
	two_auction_mechanism_generalized_destroy(TwoAuctionMechanismGeneralized *auction_mechanism)
	{
		auction_mechanism->~TwoAuctionMechanismGeneralized();
	}
}

#endif // TWO_AUCTION_GENERILIZED_WRAPPER_H
