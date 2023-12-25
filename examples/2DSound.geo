Include "2DMacro.geo";
Include "2DParam.pro";

/* Borders */

If (DomainType==0)
    Call "DomainSquare";
Else
    Call "DomainCircle";
EndIf

/* Holes */

clHoles = {};
cHoles = {};

For i In {1:nHoleCircle}
    x = XC~{i};
    y = YC~{i};
    r = RC~{i};
    Call "HoleCircle";
EndFor

For i In {1:nHoleLong}
    a = AL~{i}*Pi/180.0;
    x = XL~{i};
    y = YL~{i};
    r = RL~{i};
    w = WL~{i};
    Call "HoleLong";
EndFor

sBord = news; Plane Surface(sBord) = {clBord,clHoles[]};

If(MeshType==1)
    Transfinite Surface(sBord);
EndIf

/* Physical Entities */

Physical Surface("Domain",1000) = {sBord};

If(BcBord==1 && BcHole==1)
    Physical Curve("Absorbing",2000) = {lBord[],cHoles[]};
ElseIf(BcBord==0 && BcHole==0)
    Physical Curve("Reflecting",2000) = {lBord[],cHoles[]};
ElseIf(BcBord==1 && BcHole==0)
    Physical Curve("Absorbing",2000) = {lBord[]};
    Physical Curve("Reflecting",3000) = {cHoles[]};
ElseIf(BcBord==0 && BcHole==1)
    Physical Curve("Absorbing",2000) = {cHoles[]};
    Physical Curve("Reflecting",3000) = {lBord[]};
EndIf
