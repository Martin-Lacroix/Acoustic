Include "1DParam.pro";

/* Borders */

pDom[0] = newp; Point(pDom[0]) = {-xDom/2,0,0,LcDom};
pDom[1] = newp; Point(pDom[1]) = {xDom/2,0,0,LcDom};

lDom = newl; Line(lDom) = {pDom[0],pDom[1]};
clDom = newll; Curve Loop(clDom) = {lDom};


/* Physical Entities */

Physical Curve("Domain",1000) = {lDom};

If(BcDom==1)
    Physical Point("Absorbing",2000) = {pDom[]};
ElseIf(BcDom==0)
    Physical Point("Reflecting",2000) = {pDom[]};
EndIf
