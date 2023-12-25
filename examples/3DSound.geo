Include "3DMacro.geo";

DomainType = DefineNumber[0, Choices{0="Amphi",1="Cube"}, Highlight "Green", Name "1Mesh/1Domain"];
Bc = DefineNumber[0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "1Mesh/2Boundary"];
Lc = DefineNumber[1.0, Name "1Mesh/3Mesh size", Highlight "Green", Units "m", Min 1e-2, Max 10];

xA = DefineNumber[30.0, Name "2Amphi/1Length x", Highlight "Green", Units "m", Min 0, Max 100];
zA = DefineNumber[10.0, Name "2Amphi/2Height z", Highlight "Green", Units "m", Min 0, Max 100];
a = DefineNumber[45.0, Name "2Amphi/3Angle", Highlight "Green", Units "°", Min 0, Max 100];
c = DefineNumber[30.0, Name "2Amphi/4Crop", Highlight "Green", Units "%", Min 0, Max 100];
n = DefineNumber[2.0, Name "2Amphi/5Steps", Highlight "Green", Units "-", Min 2, Max 100];
s = DefineNumber[5.0, Name "2Amphi/6Start", Highlight "Green", Units "m", Min 0, Max 100];
e = DefineNumber[5.0, Name "2Amphi/7End", Highlight "Green", Units "m", Min 0, Max 100];
h = DefineNumber[5.0, Name "2Amphi/8Estrade Height", Highlight "Green", Units "m", Min 0, Max 100];

xC = DefineNumber[10.0, Name "3Cube/3Length x", Highlight "Green", Units "m", Min 0, Max 100];
yC = DefineNumber[10.0, Name "3Cube/4Width y", Highlight "Green", Units "m", Min 0, Max 100];
zC = DefineNumber[10.0, Name "3Cube/4Heignt z", Highlight "Green", Units "m", Min 0, Max 100];

If(DomainType==0)
    x = xA;
    z = zA;
    a = a*Pi/180;
    c = c/100.0;
    stepL = (x*(1-c)-s-e)*Cos(a)/(n-1);
    stepH = h/n;
    Call "DomainAmphi";
Else
    x = xC;
    y = yC;
    z = zC;
    Call "DomainCube";
EndIf

slDom = newsl; Surface Loop(slDom) = {sDom[]};
vDom = newv; Volume(vDom) = {slDom};

If(Bc==0)
Physical Surface("Reflecting",1000) = {sDom[]};
Else
    Physical Surface("Absorbing",1000) = {sDom[]};
EndIf
Physical Volume("Domain",2000) = {vDom};