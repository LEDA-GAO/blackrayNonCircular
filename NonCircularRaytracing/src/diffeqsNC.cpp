#include "def.h"

std::array<long double, 5> diffeqsNC(long double c1, long double c2, std::array<long double, 5> vars)
{

	long double r, th;
	long double ku, kphi;
	long double denom;
	long double t31,t32,t33,t47,t48,t49,t51,t50,t52,t53,t54,t55,t56,t57,t58,t60,t61,t62,t74,t75,t92,t90,t72,t73,t76,t79,t80,t81,t82,t83,t86, t87,t88,t89,t91;
    long double t93,t94,t96,t97,t98,t100,t99,t117,t118,t119,t120,t121,t122,t123,t124,t125,t131,t132,t149,t113,t116,t129,t130,t147,t112;   
	
    long double g_uu, g_pp, g_up, g_ur, g_rp;
    long double gurr, guur, gurp,guthth;
    long double dguudr, dgppdr, dgupdr, dgththdr, dguudth, dgppdth, dgupdth, dgrpdth, dgththdth;
    long double ch_ruu, ch_rpp, ch_rur, ch_ruth, ch_rup, ch_rrth, ch_rrp, ch_rthp,ch_rthth;
    long double ch_thuu, ch_ththth, ch_thpp, ch_thup, ch_thrth, ch_thrp;
    long double ku2, kp2, kur, kuth, kup, krth, krp, kthp, kth2;
    std::array<long double, 5> diffs;
    const int max_prec = std::numeric_limits<long double>::max_digits10;
  

	r = vars[0];
	th = vars[1];	
    t31 = pow(r,2);
    t32 = pow(spin,2);
    t33 = cos(th);
    t47 = pow(t33,2);
    t48 = t32*t47;
    t49 = t31 + t48;
    t51 = betaNC/2.;
    t50 = 1/t49;
    t52 = pow(3,t51);
    t53 = pow(4,betaNC);
    t54 = pow(lNP,4);
    t55 = pow(t49,-3);
    t56 = t54*t55;
    t57 = pow(t56,t51);
    t58 = t52*t53*t57;
    t60 = 1 + t58;
    t61 = 1/t60;
    t73 = sin(th);
    t74 = pow(t73,2);
    t91 = pow(t73,4);
    t89 = pow(spin,4);
    t62 = 2*r*t50*t61;
    t72 = -1 + t62;
    t75 = t31 + t32;
    t76 = pow(t75,2);
    t79 = -2*r*t61;
    t80 = t31 + t32 + t79;
    t81 = -(t32*t74*t80);
    t82 = t76 + t81;
    t86 = pow(r,4);
    t87 = -(t50*t74*t86);
    t88 = -2*t31*t32*t50*t74;
    t90 = -(t50*t74*t89);
    t92 = t32*t91;
    t93 = t31*t32*t50*t91;
    t94 = t50*t89*t91;
    t96 = t87 + t88 + t90 + t92 + t93 + t94;
    t97 = 1/t96;
    t99 = pow(t60,-2);
    t98 = pow(t49,-2);
    t116 = 2*betaNC;
    t117 = 1 + t116;
    t118 = pow(2,t117);
    t119 = 1 + t51;
    t120 = pow(3,t119);
    t121 = -1 + t51;
    t122 = pow(t56,t121);
    t123 = pow(t49,-5);
    t129 = 2*r;
    t130 = pow(t49,-4);
    t147 = pow(t73,3);
    t112 = pow(spin,3);
    g_ur = 1;
    g_uu = t72;
    g_pp = t50*t74*t82;
    g_rp = -(spin*t74);
    g_up = -2*r*spin*t50*t61*t74;

    gurr = t97*(t50*t72*t74*t82 - 4*t31*t32*t91*t98*t99);
    guur = -((2*t31*t32*t50*t74 + t50*t74*t86 + t50*t74*t89 - t31*t32*t50*t91 - t50*t89*t91)*t97);
    gurp = t50*(-(spin*t31*t74) - t112*t47*t74)*t97;
    guthth = t50;
    dguudr = 2*t50*t61 - 4*t31*t61*t98 + betaNC*t118*t120*t122*t123*t31*t54*t99;
    dgppdr = -2*r*t74*t82*t98 + t50*t74*(4*r*t75 - t32*t74*(t129 - 2*t61 - betaNC*t118*t120*t122*t130*t31*t54*t99));
    dgupdr = -2*spin*t50*t61*t74 + 4*spin*t31*t61*t74*t98 - spin*betaNC*t118*t120*t122*t123*t31*t54*t74*t99;
    dgththdr = t129;
    dguudth = 4*r*t32*t33*t61*t73*t98 - r*betaNC*t118*t120*t122*t123*t32*t33*t54*t73*t99;
    dgppdth = 2*t33*t50*t73*t82 + 2*t147*t32*t33*t82*t98 + t50*t74*(-2*t32*t33*t73*t80 - r*betaNC*t118*t120*t122*t130*t147*t33*t54*t89*t99);
    dgupdth = -4*r*spin*t33*t50*t61*t73 - 4*r*t112*t147*t33*t61*t98 + r*betaNC*t112*t118*t120*t122*t123*t147*t33*t54*t99;
    dgrpdth = -2*spin*t33*t73;
    dgththdth = -2*t32*t33*t73;

    ch_ruu = -0.5*gurr*dguudr;
    //std::cout <<"ch_ruu "<< std::setprecision(max_prec) << ch_ruu<< " "<<std::endl;
    ch_rthth = -0.5*gurr*dgththdr;
    //std::cout <<"ch_rthth "<< std::setprecision(max_prec) << ch_rthth<< " "<<std::endl;
    ch_rpp = -0.5*gurr*dgppdr;
    // std::cout <<"ch_rpp "<< std::setprecision(max_prec) << ch_rpp<< " "<<std::endl;
    ch_rur = 0.5*(guur*dguudr+gurp*dgupdr);
    // std::cout <<"ch_rur "<< std::setprecision(max_prec) << ch_rur<< " "<<std::endl;
    ch_ruth = 0.5*(guur*dguudth+gurp*dgupdth);
    // std::cout <<"ch_ruth "<< std::setprecision(max_prec) << ch_ruth<< " "<<std::endl;
    ch_rup = -0.5*gurr*dgupdr;
    // std::cout <<"ch_rup  "<< std::setprecision(max_prec) << ch_rup << " "<<std::endl;
    ch_rrth = 0.5*gurp*dgrpdth;
    // std::cout <<"ch_rrth "<< std::setprecision(max_prec) << ch_rrth<< " "<<std::endl;
    ch_rrp = 0.5*(guur*dgupdr+gurp*dgppdr);
    // std::cout <<"ch_rrp "<< std::setprecision(max_prec) << ch_rrp<< " "<<std::endl;
     //for ch_rthp a lot of digits are canceled out so its accuracy is not good
     //It is okay for small value. Keep in mind if future problems appear
     //possible solution: use its expression directly 
    ch_rthp = 0.5*(gurr*dgrpdth+guur*dgupdth+gurp*dgppdth);

    ch_thuu = -0.5*guthth*dguudth;
    ch_ththth = 0.5*guthth*dgththdth;
    ch_thpp = -0.5*guthth*dgppdth;
    ch_thup = -0.5*guthth*dgupdth;
    ch_thrth = 0.5*guthth*dgththdr;
    ch_thrp = -0.5*guthth*dgrpdth;


    //The part below is correct now. 
    denom = g_uu*g_pp-g_up*g_up;
	
    ku = (c1*g_pp - c2*g_up-(g_ur*g_pp-g_up*g_rp)*vars[3])/denom;
    //cout<<"ku "<<ku<<endl;
	kphi = -(c1*g_up - c2*g_uu - (g_ur*g_up-g_uu*g_rp)*vars[3])/denom;
    //cout<<"kphi "<<kphi<<endl;

    

	diffs[0] = vars[3];
	diffs[1] = vars[4];
	diffs[2] = kphi;
	
    ku2 = ku*ku;
    kp2 = kphi*kphi;
    kur = ku*vars[3];
    kuth = ku*vars[4];
    kup = ku*kphi;
    krth = vars[3]*vars[4];
    krp = vars[3]*kphi;
    kthp = vars[4]*kphi;
    kth2 = vars[4]*vars[4];

    //The second derivative for r and theta should be the same for these two coordinate
    diffs[3] = -(ch_ruu*ku2+ch_rpp*kp2+ch_rthth*kth2+2.0*(ch_rur*kur+ch_ruth*kuth+ch_rup*kup+ch_rrth*krth+ch_rrp*krp+ch_rthp*kthp));
    diffs[4] = -(ch_thuu*ku2+ch_ththth*kth2+ch_thpp*kp2+2.0*(ch_thup*kup+ch_thrth*krth+ch_thrp*krp));

    return diffs; 
}
