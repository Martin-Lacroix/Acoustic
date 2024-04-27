Macro "DomainAmphi"

  /* Stair Points */

  pStair[0] = newp; Point(pStair[0]) = {0,0,0};

  j = 1;
  k = 0;
  For i In {0:n-1}
      pStair[j+0] = newp; Point(pStair[j+0]) = {(c*x+s+i*stepL)*Cos(a),-(x*c+s+i*stepL)*Sin(a),i*stepH,Lc};
      pStair[j+1] = newp; Point(pStair[j+1]) = {(c*x+s+i*stepL)*Cos(a),(x*c+s+i*stepL)*Sin(a),i*stepH,Lc};
      pStair[j+2] = newp; Point(pStair[j+2]) = {(c*x+s+i*stepL)*Cos(a),(x*c+s+i*stepL)*Sin(a),(i+1)*stepH,Lc};
      pStair[j+3] = newp; Point(pStair[j+3]) = {(c*x+s+i*stepL)*Cos(a),-(x*c+s+i*stepL)*Sin(a),(i+1)*stepH,Lc};
      pStair[j+4] = newp; Point(pStair[j+4]) = {0,0,(i+1)*stepH};

      lStair[k+0] = newl; Circle(lStair[k+0]) = {pStair[j],pStair[j-1],pStair[j+1]};
      lStair[k+1] = newl; Line(lStair[k+1]) = {pStair[j+1],pStair[j+2]};
      lStair[k+2] = newl; Circle(lStair[k+2]) = {pStair[j+2],pStair[j+4],pStair[j+3]};
      lStair[k+3] = newl; Line(lStair[k+3]) = {pStair[j+3],pStair[j]};

      j += 5;
      k += 4;
  EndFor

  l = 0;
  For i In {0:n-2}
      lStair[k+0] = newl; Line(lStair[k+0]) = {pStair[l+4],pStair[l+6]};
      lStair[k+1] = newl; Line(lStair[k+1]) = {pStair[l+3],pStair[l+7]};

      k += 2;
      l += 5;
  EndFor

  /* Border Points */

  pBord[0] = newp; Point(pBord[0]) = {0,0,z,Lc};
  pBord[1] = newp; Point(pBord[1]) = {x*c*Cos(a),-x*c*Sin(a),0,Lc};
  pBord[2] = newp; Point(pBord[2]) = {x*Cos(a),-x*Sin(a),h,Lc};
  pBord[3] = newp; Point(pBord[3]) = {x*Cos(a),-x*Sin(a),z,Lc};
  pBord[4] = newp; Point(pBord[4]) = {x*c*Cos(a),-x*c*Sin(a),z,Lc};
  pBord[5] = newp; Point(pBord[5]) = {x*c*Cos(a),x*c*Sin(a),0,Lc};
  pBord[6] = newp; Point(pBord[6]) = {x*Cos(a),x*Sin(a),h,Lc};
  pBord[7] = newp; Point(pBord[7]) = {x*Cos(a),x*Sin(a),z,Lc};
  pBord[8] = newp; Point(pBord[8]) = {x*c*Cos(a),x*c*Sin(a),z,Lc};

  /* Amphi Lines */

  lBord[0] = newl; Line(lBord[0]) = {pBord[1],pStair[1]};
  lBord[1] = newl; Line(lBord[1]) = {pStair[j-2],pBord[2]};
  lBord[2] = newl; Line(lBord[2]) = {pBord[2],pBord[3]};
  lBord[3] = newl; Line(lBord[3]) = {pBord[3],pBord[4]};
  lBord[4] = newl; Line(lBord[4]) = {pBord[4],pBord[1]};
  lBord[5] = newl; Line(lBord[5]) = {pBord[5],pStair[2]};
  lBord[6] = newl; Line(lBord[6]) = {pStair[j-3],pBord[6]};
  lBord[7] = newl; Line(lBord[7]) = {pBord[6],pBord[7]};
  lBord[8] = newl; Line(lBord[8]) = {pBord[7],pBord[8]};
  lBord[9] = newl; Line(lBord[9]) = {pBord[8],pBord[5]};
  lBord[10] = newl; Line(lBord[10]) = {pBord[1],pBord[5]};
  lBord[11] = newl; Line(lBord[11]) = {pBord[4],pBord[8]};
  lBord[12] = newl; Circle(lBord[12]) = {pBord[2],pStair[j-1],pBord[6]};
  lBord[13] = newl; Circle(lBord[13]) = {pBord[3],pBord[0],pBord[7]};

  /* Amphi Curve Loops */

  l = 0;
  For i In {0:n-1}
      clDom[i] = newll; Curve Loop(clDom[i]) = {lStair[l],lStair[l+1],lStair[l+2],lStair[l+3]};
      l += 4;
  EndFor

  l = 0;
  For i In {0:n-2}
      clDom[i+n] = newll; Curve Loop(clDom[i+n]) = {lStair[2+l],lStair[n*4+2*i],lStair[4+l],-lStair[n*4+2*i+1]};
      l += 4;
  EndFor

  clDom[2*n-1] = newll; Curve Loop(clDom[2*n-1]) = {lBord[0],lStair[0],-lBord[5],-lBord[10]};
  clDom[2*n+0] = newll; Curve Loop(clDom[2*n+0]) = {lBord[1],lBord[12],-lBord[6],lStair[(n-1)*4+2]};
  clDom[2*n+1] = newll; Curve Loop(clDom[2*n+1]) = {lBord[4],lBord[10],-lBord[9],-lBord[11]};
  clDom[2*n+2] = newll; Curve Loop(clDom[2*n+2]) = {lBord[3],lBord[11],-lBord[8],-lBord[13]};
  clDom[2*n+3] = newll; Curve Loop(clDom[2*n+3]) = {lBord[12],lBord[7],-lBord[13],-lBord[2]};

  Stair[] = {};
  For i In {0:n-2}
      Stair[] += {-lStair[4*i+3],lStair[n*4+2*i]};
  EndFor
  Stair[] += {-lStair[4*(n-1)+3]};
  clDom[2*n+4] = newll; Curve Loop(clDom[2*n+4]) = {lBord[0],Stair[],lBord[1],lBord[2],lBord[3],lBord[4]};

  Stair[] = {};
  For i In {0:n-2}
      Stair[] += {lStair[4*i+1],lStair[n*4+2*i+1]};
  EndFor
  Stair[] += {lStair[4*(n-1)+1]};
  clDom[2*n+5] = newll; Curve Loop(clDom[2*n+5]) = {lBord[5],Stair[],lBord[6],lBord[7],lBord[8],lBord[9]};

  /* Amphi Surfaces */

  For i In{0:2*n+3}
      sDom[i] = news; Surface(sDom[i]) = {clDom[i]};
  EndFor

  sDom[2*n+4] = news; Plane Surface(sDom[2*n+4]) = {clDom[2*n+4]};
  sDom[2*n+5] = news; Plane Surface(sDom[2*n+5]) = {clDom[2*n+5]};
Return


Macro "DomainCube"

  /* Cube Points */

  pDom[0] = newp; Point(pDom[0]) = {-x/2,-y/2,-z/2,Lc};
  pDom[1] = newp; Point(pDom[1]) = {x/2,-y/2,-z/2,Lc};
  pDom[2] = newp; Point(pDom[2]) = {x/2,y/2,-z/2,Lc};
  pDom[3] = newp; Point(pDom[3]) = {-x/2,y/2,-z/2,Lc};
  pDom[4] = pDom[0];

  pDom[5] = newp; Point(pDom[5]) = {-x/2,-y/2,z/2,Lc};
  pDom[6] = newp; Point(pDom[6]) = {x/2,-y/2,z/2,Lc};
  pDom[7] = newp; Point(pDom[7]) = {x/2,y/2,z/2,Lc};
  pDom[8] = newp; Point(pDom[8]) = {-x/2,y/2,z/2,Lc};
  pDom[9] = pDom[5];

  /* Cube Lines */

  For i In{0:3}
      lDom[i] = newl; Line(lDom[i]) = {pDom[i],pDom[i+1]};
      lDom[i+5] = newl; Line(lDom[i+5]) = {pDom[i+5],pDom[i+6]};
      lDom[i+10] = newl; Line(lDom[i+10]) = {pDom[i],pDom[i+5]};
  EndFor

  /* Cube Curve Loops */

  clDom[0] = newll; Curve Loop(clDom[0]) = {lDom[0],lDom[1],lDom[2],lDom[3]};
  clDom[1] = newll; Curve Loop(clDom[1]) = {lDom[5],lDom[6],lDom[7],lDom[8]};
  clDom[2] = newll; Curve Loop(clDom[2]) = {lDom[10],lDom[5],-lDom[11],-lDom[0]};
  clDom[3] = newll; Curve Loop(clDom[3]) = {lDom[11],lDom[6],-lDom[12],-lDom[1]};
  clDom[4] = newll; Curve Loop(clDom[4]) = {lDom[12],lDom[7],-lDom[13],-lDom[2]};
  clDom[5] = newll; Curve Loop(clDom[5]) = {lDom[13],lDom[8],-lDom[10],-lDom[3]};

  /* Cube Surfaces */

  For i In{0:5}
      sDom[i] = news; Plane Surface(sDom[i]) = {clDom[i]};
  EndFor
Return
