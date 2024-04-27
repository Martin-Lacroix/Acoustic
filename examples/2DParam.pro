DefineConstant[BcBord = {0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "1Domain/1BC"}];
DefineConstant[MeshType = {0, Choices{0="Random",1="Transfinite"}, Highlight "Green", Name "1Domain/2Element"}];
DefineConstant[DomainType = {0, Choices{0="Square",1="Circle"}, Highlight "Green", Name "1Domain/3Type"}];
DefineConstant[LcBord = {1.0, Highlight "Green", Name "1Domain/4Mesh size", Units "m",Min 1e-2, Max 10}];
DefineConstant[xBord = {30.0, Highlight "Green", Name "1Domain/5Width", Units "m",Min 0, Max 100}];
DefineConstant[yBord = {30.0, Highlight "Green", Name "1Domain/6Height", Units "m", Min 0, Max 100}];
DefineConstant[rBord = {20.0, Highlight "Green", Name "1Domain/7Radius", Units "m", Min 0, Max 100}];

DefineConstant[BcHole = {0, Choices{0="Reflecting",1="Absorbing"}, Highlight "Green", Name "2Holes/1BC"}];
DefineConstant[LcHole = {1.0, Highlight "Green", Name "2Holes/2Mesh size", Units "m",Min 1e-2, Max 100}];
DefineConstant[nHoleCircle = {0, Step 1, Highlight "Green", Name "2Holes/3Circle", Units "-", Min 0, Max 20}];
DefineConstant[nHoleLong = {0, Step 1, Highlight "Green", Name "2Holes/3Long", Units "-", Min 0, Max 20}];

For i In {1:nHoleCircle}
  DefineConstant[
    XC~{i} = {0, Highlight "Green", Name Sprintf("2Holes/Circle %g/1x", i), Units "m", Min 0, Max 100},
    YC~{i} = {0, Highlight "Green", Name Sprintf("2Holes/Circle %g/2y", i), Units "m", Min 0, Max 100},
    RC~{i} = {1.0, Highlight "Green", Name Sprintf("2Holes/Circle %g/4Radius", i), Units "m", Min 0, Max 100}
  ];
EndFor

For i In {1:nHoleLong}
  DefineConstant[
    XL~{i} = {0, Highlight "Green", Name Sprintf("2Holes/Long %g/1x", i), Units "m", Min 0, Max 100},
    YL~{i} = {0, Highlight "Green", Name Sprintf("2Holes/Long %g/2y", i), Units "m", Min 0, Max 100},
    RL~{i} = {1.0, Highlight "Green", Name Sprintf("2Holes/Long %g/3Radius", i), Units "m", Min 0, Max 100},
    WL~{i} = {1.0, Highlight "Green", Name Sprintf("2Holes/Long %g/4Width", i), Units "m", Min 0, Max 100}
    AL~{i} = {0, Highlight "Green", Name Sprintf("2Holes/Long %g/5Angle", i), Units "°", Min -180, Max 180}
  ];
EndFor
