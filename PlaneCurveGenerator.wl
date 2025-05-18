(*
PlaneCurveGenerator Package
by Xah Lee
Created 1995.
Last Modified: 2024-07-17
*)

(* :Context: PlaneCurveGenerator` *)

(* :Title: PlaneCurveGenerator *)

(* :Author: Xah Lee *)

(* : Copyright 1995, 2024, Xah Lee. *)

(* :Summary:
This package exports a set of graphical functions that generate animation
showing various ways to trace special plane curves such as
conic sections, cycloids, hypotrochoids, epitrochoids , various spirals, Conchoids,
Cissoids, Witch of Angnesi... etc.
*)

(* :Keywords: graphics, curve, plot, geometry, calculus,
	normal, parallel, caustics, evolutes, radial, concoid, inversion
*)

(* :Homepage: http://xahlee.info/M/plane_curve_generator.html *)

(* :Mathematica Version: 10 *)
(* :Package Version: 2.1.20250221220942 *)
(* :History: *)

(* :Discussion:
This package exports a set of graphical functions that generate animation
showing various ways to trace special plane curves. More than 60 such functions are exported. A typical function looks like

EllipseGenerator[NumberOfFrames->30]

This package is designed to visually demonstrate various properties of special plane curves. 

 Xah
 xah@xahlee.org
 http://xahlee.org/

*)

(* :Sources:
	Mathematica Programming:
	Maeder, Roman: Programming In Mathematica, 2nd ed. chapter 1-3.
	Smith, Cameron; Blachman, Nancy: The Mathematica Graphics Guidebook.
	Mathematics of plane curves: (in order of my preference)
	Yates, R. C.; Curves and their Properties
	(The National Council of Teachers of Mathematics, 1952).
	Lockwood, E.H.: A book of Curves. (Cambridge Uni. Press, 1961).
	Lawrence, J. D.; A Catalog of Special Plane Curves (Dover, 1972).
	Seggern, D. v.; CRC Standard Curves and Surfaces (CRC Press; 1994).
	Encyclopedia Britannica, Geometry: Special Curves.
*)

(* HHHH--------------------------------------------------- *)

BeginPackage["PlaneCurveGenerator`"]

Clear[
AstroidGenerator,
AstroidGenerator2,
AstroidGenerator3,
CardioidGenerator,
ConchoidGenerator,
ConchoidOfNicomedesGenerator,
DeltoidGenerator,
DeltoidGenerator2,
DerivativeGenerator,
EllipseGenerator,
EllipseGenerator2,
EllipseGenerator3,
EllipseGenerator4,
EpitrochoidGenerator,
HyperbolaGenerator,
HypotrochoidGenerator,
LimaconOfPascalGenerator,
LimaconOfPascalGenerator2,
NephroidGenerator,
NumberOfFrames,
ParaPolarGenerator,
ParabolaGenerator,
QuadratrixOfHippiasGenerator,
RoseGenerator,
RoseGenerator2,
SineGenerator,
TractrixGenerator,
WitchOfAgnesiGenerator
];

NumberOfFrames::usage = "NumberOfFrames is an option for sevenal animation functions in the package PlaneCurveGenerator.m.
It specify the number of frames to be generated.";

ParaPolarGenerator::usage = "ParaPolarGenerator[{dist[t],angle[t]}, {t, tMin, tMax}] generates an animation tracing the curve {dist[t],angle[t]} that is parametrized for polar coordinates.
dist[t] is a function specifying it's distance to the origin at time t, angle[t] is the function specifing the angle.
Special options are NumberOfFrames->72.
Other options are passed to Graphics[].
Example: ParaPolarGenerator[{Sin[t],Cos[t]}, {t,0,2 Pi}]";

ConchoidGenerator::usage = "ConchoidGenerator[{xf[t],yf[t]}, {t,tMin,tMax,(tStep)}, {{p1,p2}, n}] generates an animation tracing the conchoid of the curve {xf[t],yf[t]} with respect to the conchoid point {p1,p2} and distance n.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
Example: ConchoidGenerator[ {t, Sin[t]}, {t,0, 3 Pi}, {{3,-3}, 2}];0;";

DerivativeGenerator::usage = "DerivativeGenerator[ f[t], {t,tMin,tMax, (tStep)}] generates an animation that trace out the derivative of f[t].
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
Example: DerivativeGenerator[Sin[x],{x,0,2 Pi}]";

AstroidGenerator::usage = "AstroidGenerator[] generates an animation of astroid by Trammel of Archimedes.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
AstroidGenerator[{tMin,tMax, (tStep)}] can also be used.
tMin >= 0 and tMax =< 2 Pi.
See ?RoseGenerator for detail.";

AstroidGenerator2::usage = "AstroidGenerator2[] generates an animation that trace out the curve astroid as a special case of hypotrochoid.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
AstroidGenerator2[{tMin,tMax, (tStep)}] can also be used.
0 <= tMin <= tMax =< 2 Pi.
See ?RoseGenerator for detail.";

AstroidGenerator3::usage = "AstroidGenerator3[] generates an animation that trace out the curve astroid showing that it is the envelope of co-axial ellipses whoes sum of major and minor axes is contsant.
Special options are NumberOfFrames->10.
Other options are passed to Graphics[].
AstroidGenerator3[{tMin, tMax, (tStep)}] can also be used.
0 <= tMin <= tMax <=1.
See ?RoseGenerator for more detail.";

CardioidGenerator::usage = "CardioidGenerator[] generates an animation that trace out the curve cardioid as a special case of epitrochoid.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
CardioidGenerator[{tMin, tMax, (tStep)}] can also be used.
0 <= tMin <= tMax <= 2 Pi.
See ?RoseGenerator for more detail.";

ConchoidOfNicomedesGenerator::usage = "ConchoidOfNicomedesGenerator[ distance1,distance2] generates an animation tracing the curve conchoid of Nicomedes (i.e.
conchoid of a line).
distance1 is the distance from the conchoid point to the line.
distance2 is the distances of the tracing points to the line.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
ConchoidOfNicomedesGenerator[d1,d2,{tMin, tMax, (tStep)}] can also be used. -Infinity < tMin < tMax < Infinity.
See ?RoseGenerator for more detail.
Example: ConchoidOfNicomedesGenerator[1,2]";

DeltoidGenerator::usage = "DeltoidGenerator[] generates an animation that trace out the curve deltoid as a special case of hypotrochoid.
Special options are NumberOfFrames->46.
Other options are passed to Graphics[].
DeltoidGenerator[{tMin, tMax, (tStep)}] can also be used.
0 <= tMin <= tMax <= 2 Pi.
See ?RoseGenerator for more detail.";

DeltoidGenerator2::usage = "DeltoidGenerator2[({t1,t2,t3})] generates an animation that trace out the curve deltoid by the envelope of Simmson lines of a triangle. {t1,t2,t3} is optional.
They specify the shape of triangle inscribed in a circle.
0 < ti < 2 Pi.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
DeltoidGenerator2[{t1,t2,t3},{tMin,tMax, (tStep)}] can also be used.
tMin >= 0 and tMax =< 2 Pi.
See ?RoseGenerator for detail.";

EllipseGenerator::usage = "EllipseGenerator[n] generates an animation that trace out an ellipse by way of glissette of one line of constant length on other two mutually orthogonal lines.
n must be a number between 0 and 1 inclusive.
The closer n is to 1/2, the less essentricity of the ellipse generated.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
EllipseGenerator[n, {tMin,tMax,(tStep)}] can also be used.
tMin and tMax control whether the animation will trace out a complete ellipse.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
tMin >=0 and tMax<= 2 Pi.
Example: EllipseGenerator[.75]";

EllipseGenerator2::usage = "EllipseGenerator2[n] generates an animation that trace out an ellipse as a special case of hypotrochoid.
n is the distance from the tracing point to the center of the rolling circle.
n>=0.
The default n is 0.2.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
EllipseGenerator2[n, {tMin,tMax,(tStep)}] can also be used.
tMin and tMax control whether the animation will trace out a complete ellipse.
tMin >=0 and tMax<= 2 Pi.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
Example: EllipseGenerator2[.2 (*try 0,.5, and 2 *)]";

EllipseGenerator3::usage = "EllipseGenerator3[e] generates an animation that trace out the curve ellipse similating a string of contant length tied to two tacks.
e is the essentricity of the ellipse.
0 < e < 1.
The default of e is .8.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
EllipseGenerator3[e, {tMin, tMax, (tStep)}] can also be used.
0 <= tMin <= tMax <= 2 Pi.
See ?RoseGenerator for more detail.
Example: EllipseGenerator3[.8]";

EllipseGenerator4::usage = "EllipseGenerator4[{a, b}] generates an animation that trace out the curve ellipse by a mechanism implicit in the parametrization {a Cos[t], b Sin[t]}.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
EllipseGenerator4[{a,b},{tMin, tMax, (tStep)}] can also be used.
0 <= tMin < tMax =< 2 Pi.
See ?RoseGenerator for detail.
To get ellipse with eccentricity e, let b == Sqrt[a - a e^2].
Example: EllipseGenerator4[{2,1},{0,4}, NumberOfFrames->50]";

EpitrochoidGenerator::usage = "EpitrochoidGenerator[{a,b,h},({tMin,tMax,(tStep)})] generates an animation that trace out an epitrochoid by means of roulette of two circles.
a and b are the radii of the fixed and rolling circle respectively.
h is the distance from the tracing point to the center of the rolling circle.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
The optional list {tMin, tMax, (tStep)} control the number of rotations, which determine whether the animation will complete tracing the curve.
0 <= tMin <= tMax < Infinity.
It is useful if you don't have enough memory to complete an animation in one session.
You can set tMin in the second session to whatever value you stopped previously.
The default is {0, n 2 Pi}, where n is an integer calculated by an internal algorithm such that the epitrochoid is closed.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
Example: EpitrochoidGenerator[ {1,2/5,1}, NumberOfFrames->260]";

HyperbolaGenerator::usage = "HyperbolaGenerator[{a, b}] generates an animation that trace out hyperbola by a mechanism implicit in the parametrization {a Sec[t], b Tan[t]}.
a > 0 and b > 0.
Default {a,b} is {2,1}.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
HyperbolaGenerator[{a,b},{tMin, tMax, (tStep)}] can also be used.
0 <= tMin < tMax <= Pi/2.
See ?RoseGenerator for detail.
Example: HyperbolaGenerator[{2,1},{0,1.4}]";

HypotrochoidGenerator::usage = "HypotrochoidGenerator[{a,b,h},({tMin,tMax,(tStep)})] generates an animation that trace out a Hypotrochoid by means of roulette of two circles.
a is the radius of the fixed circle.
b is the radius of the rolling circle.
a > b.
h is the distance from the tracing point to the center of the rolling circle.
Special options are NumberOfFrames->41.
Other options are passed to Graphics[].
The optional list {tMin, tMax, (tStep)} control the number of rotations.
If not specified, an internal algorithm will figure out the number of rotations so that the curve is closed.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
0 <= tMin <= tMax < Infinity.
Example: HypotrochoidGenerator[ {1,2/5,1}, NumberOfFrames->260]";

LimaconOfPascalGenerator::usage = "LimaconOfPascalGenerator[n] generates an animation that trace out the curve Limacon of Pascal as a special case of epitrochoid.
n is the distance from the tracing point to the rolling unit circle.
Special options are NumberOfFrames->60.
Other options are passed to Graphics[].
LimaconOfPascalGenerator[ n, {tMin, tMax,(tStep)}] can also be used.
See ?RoseGenerator for similiar usage.
Example: LimaconOfPascalGenerator[.6 (*try 0, 1, 2 *)]";

LimaconOfPascalGenerator2::usage = "LimaconOfPascalGenerator2[k] generates an animation tracing the curve limacon of Pascal as a conchoid of a unit circle with conchoid point on the circle.
k is the distances of the tracing points to the circle.
k >= 0.
Special options are NumberOfFrames->40.
Other options are passed to Graphics[].
LimaconOfPascalGenerator2[ d, {tMin, tMax,(tStep)}] can also be used.
0 <= tMin <= tMax <= 2 Pi.
See ?RoseGenerator for similiar usage.
Example: LimaconOfPascalGenerator2[ 2 (*try 0,1, or 3*)]";

NephroidGenerator::usage = "NephroidGenerator[] generates an animation that trace out the curve nephroid as a special case of epitrochoid.
Special options are NumberOfFrames->60.
Other options are passed to Graphics[].";

ParabolaGenerator::usage = "ParabolaGenerator[] generates an animation showining that the distance from points on parabola to its focus and to its directrix is equal.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
ParabolaGenerator[{tMin,tMax, (tStep)}] can also be used.
tMin and tMax specifies the starting and ending points of the parabola to be traced, they correspond to the parameter t in {t, 1/4 t^2}.
tStep effects the number of frames indirectly and overrides NumberOfFrames option.";

QuadratrixOfHippiasGenerator::usage = "QuadratrixOfHippiasGenerator[] generates an animation showining the mechanism to draw Quadratrix of Hippias.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
QuadratrixOfHippiasGenerator[{tMin,tMax, (tStep)}] can also be used.
0.01 <= tMin <= tMax <= Pi/2.
See ?RoseGenerator for similiar usage.
Example: QuadratrixOfHippiasGenerator[]";

RoseGenerator::usage = "RoseGenerator[n] generates an animation that trace out the curve rose as a special case of hypotrochoid.
n must be an integer greater than 1.
If n is even, the rose will have 2 n pedals.
If n is odd, it will have n pedals.
The default of n is 3.
Special options are NumberOfFrames->61.
Other options are passed to Graphics[].
RoseGenerator[ n, {tMin, tMax,(tStep)}] can also be used.
tMin and tMax control whether the animation will trace out a complete rose.
It is useful if you don't have enough memory to complete an animation in one session.
You can set tMin in the second session to whatever value you stopped previously.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
If n is odd, {0,Pi} will complete a rose.
If n is even, you need {0, 2 Pi}.
Example: RoseGenerator[3]";

RoseGenerator2::usage = "RoseGenerator2[n] generates an animation that trace out the curve rose by a mechanism implicit in the polar coordinate.
n must be an integer greater than 1.
If n is even, the rose will have 2 n pedals.
If n is odd, it will have n pedals.
The default of n is 3.
Special options are NumberOfFrames->72.
Other options are passed to Graphics[].
RoseGenerator2[ n, {tMin, tMax,(tStep)}] can also be used.
tMin and tMax control whether the animation will trace out a complete rose.
It is useful if you don't have enough memory to complete an animation in one session.
You can set tMin in the second session to whatever value you stopped previously.
tStep control the number of frames indirectly and overrides the NumberOfFrames option.
If n is odd, {0,Pi} will complete a rose.
If n is even, you need {0, 2 Pi}.
Example: RoseGenerator2[3]";

SineGenerator::usage = "SineGenerator[] generates an animation illustrating Sine curve as a circular function.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
SineGenerator[{tMin, tMax,(tStep)}] can also be used.
0 <= tMin < tMax <= Pi/2.
See ?RoseGenerator for detail.";

TractrixGenerator::usage = "TractrixGenerator[] generates an animation that trace out the curve Tractrix by dragging a line of fixed length.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
TractrixGenerator[{tMin,tMax, (tStep)}] can also be used.
0 <= tMin <= tMax < Pi/2.
See ?RoseGenerator for detail.
Example: TractrixGenerator[{0,1.5}]";

WitchOfAgnesiGenerator::usage = "WitchOfAgnesiGenerator[] generates an animation that trace out the curve Witch of Agnesi.
Special options are NumberOfFrames->30.
Other options are passed to Graphics[].
WitchOfAgnesiGenerator[{tMin,tMax, (tStep)}] can also be used.
tMin and tMax specifies the starting and ending points of the curve to be traced. -Pi/2 < tMin < tMax < Pi/2.
tStep effects the number of frames indirectly and overrides NumberOfFrames option.";

(* HHHH--------------------------------------------------- *)

Begin["`Private`"]

(* getPedalPoint takes three vectors cT, tp, and p. It returns a vector of the form {x,y} such that {x,y} is parallel to cT and passing tp, and perpendicular to the vector {x,y}-p. *)

Clear[ getPedalPoint ];

getPedalPoint[cT_, tp_, p_]:=
	Module[{cT2 = Normalize[cT]},
		If[ cT2==={0,0}, {}, tp + (p - tp).cT2 cT2 ]
	]

(* HHHH--------------------------------------------------- *)

Options[ParaPolarGenerator] =
	Join[
		{NumberOfFrames->72},
		Options[Graphics]
	];

ParaPolarGenerator[
	{expx_,expy_},{t_Symbol, tmin_, tmax_, dt_:Automatic},opts:OptionsPattern[]
]:=
ParaPolarGenerator[ Function@@{{t},{expx,expy}},{tmin,tmax, dt},opts];

ParaPolarGenerator[ curve_Function, {tmin_,tmax_, dt_:Automatic},
	opts:OptionsPattern[]] :=
Module[{nOfFrames,t, tMin,tMax,tStep,
	symb1, symb2, distFun, angleFun, xCoordFun, yCoordFun,
	tRangeList, pointsInRectC, pointsInPolarC, distMax,
	staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax} = N[{tmin, tmax}];
	tStep = If[ dt === Automatic, N[(tMax-tMin)/(nOfFrames-1)], N@dt ];

	tRangeList = Range[tMin, tMax, tStep];

	distFun = Compile@@{{t}, First@curve@t};
	angleFun = Compile@@{{t}, Last@curve@t};

	xCoordFun = Compile@@{{t}, First@curve@t Cos@Last@curve@t};
	yCoordFun = Compile@@{{t}, First@curve@t Sin@Last@curve@t};

	pointsInPolarC = {distFun@#, angleFun@#}& /@ tRangeList;
	pointsInRectC = {xCoordFun@#, yCoordFun@#}& /@ tRangeList ;

	distMax = Max@ Abs@ First@ Transpose@ pointsInPolarC;

	staticGP = {Hue[.65,1,.7], Thickness[.006], Circle[{0,0}, distMax]};
	movingGP =
		{Hue[.65,1,.7], Thickness[.007],
			Line[{
				{0,0},
				distMax *
				Normalize[{
					Cos@ Last@ pointsInPolarC[[#]],
					Sin@ Last@ pointsInPolarC[[#]]
				}]
			}],
		If[ Sign@ First@ pointsInPolarC[[#]] === 1,
			Hue[.3,.7,.85],
			Hue[.18,1,.92]
		],
			Thickness[.013], Line[{{0,0}, pointsInRectC[[#]] }],
		Hue[.65,1,.7], PointSize[.04], Point[{0,0}],
		Hue[0], PointSize[.02], Point@ pointsInRectC[[#]]
		}&;

	lastFrameGP = {staticGP, movingGP @ Length @ tRangeList,
		Hue[0], Point /@ pointsInRectC};

		Table[
				Graphics[
					{	staticGP, movingGP[t],
						Hue[0], Point /@ Take[ pointsInRectC,t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->All
				],
			{t, Length@tRangeList}
		]
]

(* HHHH--------------------------------------------------- *)

Options[ConchoidGenerator] =
	Join[
		{NumberOfFrames->40},
		Options[Graphics]
	];

ConchoidGenerator[
	{expx_,expy_},{t_Symbol, tmin_, tmax_, dt_:Automatic},
	{{p1_,p2_}, n_}, opts:OptionsPattern[]
]:=
ConchoidGenerator[ Function@@{{t},{expx,expy}},{tmin,tmax, dt},
	{{p1,p2},n}, opts
];

ConchoidGenerator[ curve_Function, {tmin_,tmax_, dt_:Automatic},
	{{p1_,p2_}, n_},
	opts:OptionsPattern[]] :=
Module[{nOfFrames, t, tMin, tMax, tStep,
	symb1, symb2, xCoordFun, yCoordFun,
	tRangeList, curvePoints,conchoidPoints1,conchoidPoints2,
	xCoordMin,xCoordMax, yCoordMin, yCoordMax,
	pointListGP,staticGP, movingGP, lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax} = N[{tmin, tmax}];
	tStep = If[ dt === Automatic, N[(tMax-tMin)/(nOfFrames-1)], N@dt];

	tRangeList = Range[tMin, tMax, tStep];

	xCoordFun = Compile@@{{t}, First@ curve@ t};
	yCoordFun = Compile@@{{t}, Last@ curve@ t};

	curvePoints = {xCoordFun@#, yCoordFun@#}& /@ tRangeList ;
	conchoidPoints1 = (n Normalize[#-{p1,p2}])& /@ curvePoints;
	conchoidPoints2 = -conchoidPoints1;
	conchoidPoints1 = conchoidPoints1 + curvePoints;
	conchoidPoints2 = conchoidPoints2 + curvePoints;

	symb1 = Flatten[{
			First@Transpose@curvePoints,
			First@Transpose@conchoidPoints1,
			First@Transpose@conchoidPoints2
	}];

	symb2 = Flatten[{
			Last@Transpose@curvePoints,
			Last@Transpose@conchoidPoints1,
			Last@Transpose@conchoidPoints2
	}];

	{xCoordMin, xCoordMax} = {Min[symb1,p1],Max[symb1,p1]};
	{yCoordMin, yCoordMax} = {Min[symb2,p2],Max[symb2,p2]};

	staticGP = {Hue[.4,1,.5],
		First@ ParametricPlot[ Evaluate@curve@t, {t,tMin,tMax}]
	};
	movingGP =
		{ (*Hue[.18], Thickness[.004], Circle[curvePoints[[#]],n],*)
		Hue[.51,.9,.8], Thickness[.0058],
			Line[{{p1,p2}, conchoidPoints2[[#]]}],
		Hue[.6,.8,.9], Thickness[.006],
			Line[{conchoidPoints1[[#]], conchoidPoints2[[#]]}],
		Hue[.4,1,.5], PointSize[.02], Point@ curvePoints[[#]],
		Hue[.65,1,.7], PointSize[.02], Point[{p1,p2}],
		Hue[.75,1,.7], PointSize[.02], Point@ conchoidPoints1[[#]],
		Hue[0,1,1], PointSize[.02], Point@ conchoidPoints2[[#]]
		}&;

	lastFrameGP = {staticGP, movingGP@ Length@ tRangeList,
		Hue[.75,1,.7], PointSize[.008], Point /@ conchoidPoints1,
		Hue[0], PointSize[.008], Point /@ conchoidPoints2
	};

		Table[
				Graphics[
					{	staticGP, movingGP[t],
						Hue[.75,1,.7], PointSize[.008],
							Point/@ Take[conchoidPoints1, t],
						Hue[0], PointSize[.008],
							Point/@ Take[conchoidPoints2, t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->1.05 {{xCoordMin, xCoordMax},{yCoordMin,yCoordMax}}
				],
			{t, Length@tRangeList}
		]
]

(* HHHH--------------------------------------------------- *)

Options[DerivativeGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

DerivativeGenerator[ expr_,
	{t_Symbol, tmin_, tmax_, dt_:Automatic},opts:OptionsPattern[]
] :=
	DerivativeGenerator[ Function[{t},expr],
		{tmin,tmax,dt}, opts
	]

DerivativeGenerator[ fun_Function, tRange_List, opts:OptionsPattern[]] :=
Module[{nOfFrames,t,tMin,tMax,tStep,
	cP, dP, sf, tRangeList,
	staticGP, movingTangentGP, movingLinesGP, movingPointsGP,
	lastFrameGP, pr (* plot range *)},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Last@ tRange === Automatic,
			N@{First@ tRange, tRange[[2]], (tMax-tMin)/(nOfFrames-1) },
			N@ tRange
		];

	(sf = 1);(* Scaling Factor. It's there for possible future modification *)
	tRangeList = Range[tMin,tMax,tStep];
	cP = N[ fun /@ tRangeList];
	dP = N[ Derivative[1][fun] /@ tRangeList ];

	staticGP = First@ Plot[ fun[t], {t,tMin,tMax},
		PlotStyle->{{Hue[.65,1,.8],Thickness[.005]}}];

(* movingGP(s) are functions mapping into integers *)

	movingTangentGP =
	{	Hue[.35,1,.45], Thickness[.004],
		Line[{ -2 sf{1,dP[[#]]} + {tRangeList[[#]],cP[[#]]},
				2 sf{1,dP[[#]]} + {tRangeList[[#]],cP[[#]]}
		}],
		Line[{ -2 sf{1,dP[[#]]} + {tRangeList[[#]],0},
				2 sf{1,dP[[#]]} + {tRangeList[[#]],0}
		}]
	}&;

	movingLinesGP =
	{	Hue[.45,1,.45], Thickness[.004],
		Line[{
			{ tRangeList[[#]],
				If[	Sign@ cP[[#]] === 1,
					Min[ dP[[#]] sf, 0],
					Max[ dP[[#]] sf, 0]
				]
			},
			{tRangeList[[#]],
				If[	Sign@ cP[[#]] === 1,
					Max[ cP[[#]] sf, dP[[#]] sf ],
					Min[ cP[[#]] sf, dP[[#]] sf ]
				]
			}
		}],
		Line[{
			{tRangeList[[#]] + sf, dP[[#]] sf},
			{tRangeList[[#]] + sf, dP[[#]] sf +cP[[#]]}
		}],
		Hue[.45,1,.4], Thickness[.0065],
		Line[{
			{tRangeList[[#]], dP[[#]] sf},
			{tRangeList[[#]] + sf, dP[[#]] sf}
		}]
	}&;

	movingPointsGP =
	{	PointSize[.02], Hue[.35,1,.5],
		Point[{tRangeList[[#]], cP[[#]]}],
		Point[{tRangeList[[#]]+sf, dP[[#]] sf }],
		Hue[0,1,.9], Point[{tRangeList[[#]], 0}],
		Hue[0,1,.9], Point[{tRangeList[[#]], dP[[#]] sf }]
	}&;

	lastFrameGP = {staticGP,
		Hue[0], PointSize[.01],
		Point /@ Transpose[{tRangeList, dP}],
		movingTangentGP@ Length@ tRangeList,
		movingLinesGP@ Length@ tRangeList,
		movingPointsGP@ Length@ tRangeList
	};

	pr = If[
		MemberQ[{cP,dP},-Infinity|Infinity|ComplexInfinity|Indeterminate,-1],
		{{tMin, tMax+sf} 1.05, {tMin, tMax+sf} 1.05},
		{{tMin, tMax+sf} 1.05, {Min[cP,dP],Max[cP,dP]} 1.1}
	];

		Table[
				Graphics[
					{ staticGP,
						{Hue[0,1,.9],
							Point/@ Take[ Transpose[{tRangeList, dP}],t]},
						movingTangentGP[t],
						movingLinesGP[t],
						movingPointsGP[t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				Axes->True,
				PlotRange->pr
				],
			{t, Length@ tRangeList}
		]
]

(* HHHH--------------------------------------------------- *)

Options[AstroidGenerator] =
	Join[
		{NumberOfFrames->41},
		Options[Graphics]
	];

AstroidGenerator[tRange_List:{0, 2 Pi}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,
	tRangeList, lineListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin, tMax, tStep];
	lineListGP = Line[{{Cos@#,0},{0,Sin@#}}]& /@ tRangeList;
	staticGP = {Hue[.75], Line[{{-1,0},{1,0}}],Hue[.3,1,.5],Line[{{0,-1},{0,1}}]};
	movingGP =
		{	Hue[0], Thickness[.008], Line[{{Cos@#,0},{0,Sin@#}}],
			Hue[.75], PointSize[.02], Point[{Cos@#,0}],
			Hue[.3,1,.5], Point[{0,Sin@#}]
		}&;
	lastFrameGP = {Hue[0,.5,1], lineListGP, staticGP, movingGP@tMax};

Table[
 Graphics[ {Hue[0,.5,1], Take[ lineListGP,t], staticGP, movingGP@tRangeList[[t]] },
FilterRules[ {opts}, Options[ Graphics ] ],
AspectRatio->Automatic, Axes->True,
Ticks->{{-1,-.5,.5,1},{-1,-.5,.5,1}},
PlotRange->{{-1.05,1.05},{-1.05,1.05}}
],
			{t, Length@tRangeList}
		]
]

AstroidGenerator2[tRange_List:{Automatic}, opts:OptionsPattern[]] := HypotrochoidGenerator[{1,1/4, 1/4}, tRange, opts, NumberOfFrames -> 81]

Options[AstroidGenerator3] =
	Join[
		{NumberOfFrames->10},
		Options[Graphics]
	];

AstroidGenerator3[tRange_List:{0, 1}, opts:OptionsPattern[]] :=
Module[{nOfFrames,a,b,i,t,tMin,tMax,tStep,
	tRangeList, ellipseListGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange === 3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin, tMax, tStep];
	ellipseListGP = Table[
		First@ParametricPlot[ {i Cos@t, (1-i) Sin@t}, {t,0, 2 Pi}],
		{i, tMin, tMax, tStep}
	];

	lastFrameGP = {Hue[0,.5,1], ellipseListGP,
		Hue[0], Thickness[.005], Last@ellipseListGP
	};

		Table[
				Graphics[
					{Hue[0,.5,1], Take[ ellipseListGP,t],
					Hue[0], Thickness[.005], ellipseListGP[[t]]
					},
FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				Ticks->{{-1,-.5,.5,1},{-1,-.5,.5,1}},
				PlotRange->1.05 {{-1,1},{-1,1}}
				],
			{t, Length@tRangeList}
		]
]

CardioidGenerator[ tRange_List:{Automatic}, opts:OptionsPattern[]] :=
EpitrochoidGenerator[{1,1,1}, tRange, opts, NumberOfFrames->40]

ConchoidOfNicomedesGenerator[
	dist1_:1, dist2_:2, tRange_List:{Automatic}, opts:OptionsPattern[]
] :=
	If[ tRange==={Automatic},
		ConchoidGenerator[ {#,dist1}&,(dist1+dist2) {-1,1},{{0,0},dist2},opts],
		ConchoidGenerator[ {#,dist1}&,tRange,{{0,0},dist2},opts]
	]/;NumberQ@ N@ dist1 && NumberQ@ dist2

DeltoidGenerator[tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	HypotrochoidGenerator[{1,1/3, 1/3}, tRange, opts, NumberOfFrames->46]

Options[DeltoidGenerator2] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

DeltoidGenerator2[
	triangle_List:{Automatic},
	tRange_List:{0, 2 Pi}, opts:OptionsPattern[]
] :=
Module[{pA,pB,pC,p1,p2,p3,
	nOfFrames,t,tMin,tMax,tStep,
	tRangeList, movingAndTrailListGP,trailListGP,
	staticGP, lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@ tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	{pA,pB,pC} = If[ triangle==={Automatic},
		N@ {Cos@#, Sin@#}& /@ {1,1.7,3.5},
		N@ {Cos@#, Sin@#}& /@ triangle
	];

	tRangeList = Range[tMin, tMax, tStep];

	staticGP =
	{	Hue[0,0,0], Circle[{0,0},1],
		Thickness[.006], Line[{pA,pB,pC,pA}],
		GrayLevel[.5], PointSize[.03], Point /@ {pA,pB,pC}
	};

	movingAndTrailListGP =
	{	p1 = getPedalPoint[ pB-pA, pA, #];
		p2 = getPedalPoint[ pC-pB, pB, #];
		p3 = getPedalPoint[ pA-pC, pC, #];
		Thickness[.004], Hue[0,0,0],
		Line[{ #, p1, pB, pA}],
		Line[{ #, p2, pC, pB}],
		Line[{ #, p3, pA, pC}],
		{Hue[.65], PointSize[.02], Point /@ {p1,p2,p3}  },
		{Hue[0], PointSize[.02], Point[#]},
		(* Simson line *)
		Hue[0], Thickness[.007],
		Line[
			If[ Sqrt[(p1-p2).(p1-p2)] < .1,
				{(p1-p3) 50 + p3, (p1-p3) (-50) + p3},
				{(p1-p2) 50 + p3, (p1-p2) (-50) + p3}
			]
		]
	}& /@ ( {Cos@#, Sin@#}& /@ tRangeList );

	(*a list of Simson lines *)
	trailListGP = Last/@ movingAndTrailListGP;

	lastFrameGP = {
		Hue[0,.5,1], Drop[ trailListGP, -1],
		Last@ movingAndTrailListGP, staticGP
	};

		Table[
				Graphics[
					{	Hue[0,.5,1], Take[ trailListGP, t-1],
						movingAndTrailListGP[[t]],
						staticGP
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic,
				PlotRange->{{-1.05,1.05},{-1.05,1.05}} 2
				],
			{t, Length@ tRangeList}
		]
]

Options[EllipseGenerator] =
	Join[
		{NumberOfFrames->41},
		Options[Graphics]
	];

EllipseGenerator[n_:.7, tRange_List:{0, 2 Pi}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,
	tRangeList, pointListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin, tMax, tStep];
	pointListGP =
		Point[({Cos@#,0}-{0,Sin@#}) /
				Sqrt[Plus@@(({Cos@#,0}-{0,Sin@#})^2)]*
				n + {0,Sin@#}
		]& /@ tRangeList;
	staticGP = {Hue[.75], Line[{{-1,0},{1,0}}],
		Hue[.3,1,.5],Line[{{0,-1},{0,1}}]
	};
	movingGP =
		{	Hue[.65,1,.7], Thickness[.008], Line[{{Cos@#,0},{0,Sin@#}}],
			PointSize[.02],
			Hue[.75], Point[{Cos@#,0}],
			Hue[.3,1,.5], Point[{0,Sin@#}],
			Hue[0],
			Point[({Cos@#,0}-{0,Sin@#}) /
				Sqrt[Plus@@(({Cos@#,0}-{0,Sin@#})^2)]*
				n + {0,Sin@#}]
		}&;

(* need to rewrite this so the tracing point in movingGP is done
	by the parametric ellipse formula instead by mechanical way.xxxxx
*)
	lastFrameGP = {Hue[0], pointListGP, staticGP, movingGP@tMax};

		Table[
				Graphics[
					{Hue[0], Take[ pointListGP,t],
						staticGP,
						movingGP @ tRangeList[[t]]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				Ticks->{{-1,-.5,.5,1},{-1,-.5,.5,1}},
				PlotRange->{{-1.05,1.05},{-1.05,1.05}}
				],
			{t, Length@tRangeList}
		]

]/; NumberQ@N[n]

EllipseGenerator2[n_:.2, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	HypotrochoidGenerator[{1,1/2, n}, tRange, opts] /; NumberQ[N@n]

Options[EllipseGenerator3] =
	Join[
		{NumberOfFrames->40},
		Options[Graphics]
	];

EllipseGenerator3[e_:.8, tRange_List:{0, 2 Pi}, opts:OptionsPattern[]] :=
Module[{nOfFrames,a,t,tMin,tMax,tStep,
	tRangeList, pointListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];
	tRangeList = Range[tMin, tMax, tStep];

	pointListGP = Point[{ (a=1) Cos@#, Sqrt[a - a e^2] Sin@#}]& /@ tRangeList;
	staticGP = {PointSize[.02], Hue[.75], Point[{-a e,0}], Point[{a e,0}] };
	movingGP =
		{	curvePoint = {a Cos@#, Sqrt[1-e^2] Sin@#};
			Hue[.65,1,.7], Thickness[.004],
				Line[{{-a e,0},curvePoint}],
				Line[{{ a e,0},curvePoint}],
			Hue[.0], PointSize[.02], Point[curvePoint]
		}&;

	lastFrameGP = {Hue[0], pointListGP, movingGP@tMax, staticGP};

		Table[
				Graphics[
					{Hue[0], Take[ pointListGP,t],
						staticGP,
						movingGP @ tRangeList[[t]]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				Ticks->{{-1,-.5,.5,1},{-1,-.5,.5,1}},
				PlotRange->1.05 {{-a, a},{-Sqrt[1-e^2],Sqrt[1-e^2]}}
				],
			{t, Length@tRangeList}
		]
]/; (NumberQ[N@e] && e > 0 && e < 1)

Options[EllipseGenerator4] =
	Join[
		{NumberOfFrames->40},
		Options[Graphics]
	];

EllipseGenerator4[ radii_List:{2,1}, tRange_List:{0, 2 Pi}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep, maxR, minR,
	tRangeList, pointList, trailPointList,
	staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	{minR, maxR} = Sort@ radii;
	tRangeList = Range[tMin, tMax, tStep];
	pointList = {Cos@#, Sin@#}& /@ tRangeList;
	trailPointList = {maxR, minR} # & /@ pointList;

	staticGP = {Hue[0,0,0], Circle[{0,0}, maxR], Circle[{0,0}, minR],
		PointSize[.02], Point[{0,0}]
	};
(* movingGP is a function with integer arguments *)
	movingGP =
	{	Hue[.65], Thickness[.008],
		Line[{ {0,0}, maxR pointList[[#]]}],
		Thickness[.004],
		Hue[.4, 1, .5], Line[{ maxR pointList[[#]], trailPointList[[#]] }],
		Line[{ trailPointList[[#]], minR pointList[[#]] }],
		Hue[0], PointSize[.02], Point@ trailPointList[[#]]
	}&;

	lastFrameGP = {staticGP,
		Hue[0], PointSize[.01], Point /@ trailPointList,
		movingGP@ Length@ tRangeList
	};

		Table[
				Graphics[
					{ staticGP,
					Hue[0], PointSize[.01],Point/@ Take[ trailPointList,t],
						movingGP[t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->{{-1,1},{-1,1}} maxR 1.05
				],
			{t, Length@ tRangeList}
		]
]

Options[EpitrochoidGenerator] =
	Join[
		{NumberOfFrames->41},
		Options[Graphics]
	];

EpitrochoidGenerator[{a_,b_,h_}, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep, margin,
	tRangeList, curveXValueFun, curveYValueFun,center,
	pointListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		Which[
			tRange === {Automatic},
				{	0,
					N[ Numerator@Rationalize[N@b,0] 2 Pi],
					N[(tMax-tMin)/(nOfFrames-1)]
				},
			Length@tRange === 3,
				N@tRange,
			True,
				N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1)}
		];

	margin = (a + 2 b + Max[0, h - b]) 1.05;
	tRangeList = Range[tMin, tMax, tStep];

	curveXValueFun = Compile[{t},(a+b) Cos[t] + h Cos[(a+b)/b t] ];
	curveYValueFun = Compile[{t},(a+b) Sin[t] + h Sin[(a+b)/b t] ];

	pointListGP =
		Point[{curveXValueFun@#, curveYValueFun@# }]& /@ tRangeList;

	staticGP = {Hue[.65,1,.7], Thickness[.007], Circle[{0,0}, a]};
	movingGP =
		{	center = (a+b) {Cos[tMin + (#-1) tStep], Sin[tMin + (#-1) tStep]};
			GrayLevel[.6], Disk[ center,b],
			Hue[.65,1,.7], Thickness[.008], Circle[center,b],
			Line[{ center, pointListGP[[#,1]] }],
			Hue[0], PointSize[.02], pointListGP[[#]]
		}&;

	lastFrameGP = {staticGP, movingGP@Length@tRangeList,Hue[0],pointListGP};

		Table[
				Graphics[
					{	staticGP, movingGP[t],
						Hue[0], Take[ pointListGP,t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->{{-margin,margin},{-margin,margin}}
				],
			{t, Length@tRangeList}
		]
]

Options[HyperbolaGenerator] =
 Join[
  {NumberOfFrames->30},
  Options[Graphics]
 ];

HyperbolaGenerator[ radii_List:{2,1}, tRange_List:{0, 1.4}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep, a,b, maxR,
 tRangeList, sec, tan, trailPointList,
 staticGP, movingGP,lastFrameGP},

 nOfFrames = OptionValue[ NumberOfFrames ];

 {tMin, tMax, tStep} =
  If[ Length@tRange===3,
   N@tRange,
   N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
  ];

 {a,b} = radii; maxR = Max[ a,b];
 tRangeList = Table[ Sqrt[t], {t, tMin^2, tMax^2, (tMin+tMax) tStep} ];
 sec = Sec /@ tRangeList;
 tan = Tan /@ tRangeList;
 trailPointList = MapThread[{a #1,b #2}&,{sec,tan}];

 staticGP = {Hue[0,0,0], Circle[{0,0}, a], Circle[{0,0}, b],
  PointSize[.02], Point[{0,0}],
  Line[{{a,-50 a},{a, 50 a}}],
  Line[{{b,-50 b},{b, 50 b}}]
 };

(* movingGP is a function for mapping into integers *)
 movingGP =
 { Hue[.65], Thickness[.008],
  Line[{ {0,0}, maxR {1,1 tan[[#]]} }],
  Hue[.4, 1, .5], Thickness[.004],
  Circle[{0,0}, a sec[[#]], {0, tRangeList[[#]] }],
  Line[{
   {a sec[[#]],0},
   {a sec[[#]], b tan[[#]] },
   {b, b tan[[#]]}
  }],
  PointSize[.02],
  Point[{a, a tan[[#]]}], Point[{b,b tan[[#]]}],
  Hue[0], Point@{a sec[[#]], b tan[[#]] }
 }&;

 lastFrameGP = {staticGP,
  Hue[0], PointSize[.01], Point /@ trailPointList,
  movingGP@ Length@ tRangeList
 };

  Table[
    Graphics[
     { staticGP,
     Hue[0], PointSize[.01],Point/@ Take[ trailPointList,t],
      movingGP[t]
     }, FilterRules[ {opts}, Options[ Graphics ] ],
    AspectRatio->Automatic, Axes->True,
    PlotRange->{{-maxR, a 6},{-maxR, b 6}}
    ],
   {t, Length@ tRangeList}
  ]
]

Options[HypotrochoidGenerator] =
 Join[
  {NumberOfFrames->41},
  Options[Graphics]
 ];

HypotrochoidGenerator[{a_,b_,h_}, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep, margin,
	tRangeList, curveXValueFun, curveYValueFun,center,
	pointListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		Which[
			tRange === {Automatic},
				{	0,
					N[ Numerator@Rationalize[N@b,0] 2 Pi],
					N[(tMax-tMin)/(nOfFrames-1)]
				},
			Length@tRange === 3,
				N@tRange,
			True,
				N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1)}
		];

	margin = (a + Max[0, h - b]) 1.05;
	tRangeList = Range[tMin, tMax, tStep];

	curveXValueFun = Compile[{t},(a-b) Cos[t] + h Cos[-(a-b)/b t] ];
	curveYValueFun = Compile[{t},(a-b) Sin[t] + h Sin[-(a-b)/b t] ];

	pointListGP =
		Point[{curveXValueFun@#, curveYValueFun@# }]& /@ tRangeList;

	staticGP = {Hue[.65,1,.7], Thickness[.007], Circle[{0,0}, a]};
	movingGP =
		{	center = (a-b) {Cos[tMin + (#-1) tStep], Sin[tMin + (#-1) tStep]};
			Hue[.17], Disk[ center, b],
			Hue[.65,1,.7], Thickness[.008], Circle[center,b],
			Line[{ center, pointListGP[[#,1]] }],
			Hue[0], PointSize[.02], pointListGP[[#]]
		}&;

	lastFrameGP = {staticGP, movingGP@Length@tRangeList,Hue[0],pointListGP};

		Table[
				Graphics[
					{	staticGP, movingGP[t],
						Hue[0], Take[ pointListGP,t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->{{-margin,margin},{-margin,margin}}
				],
			{t, Length@tRangeList}
		]
]

LimaconOfPascalGenerator[ n_:.6, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	EpitrochoidGenerator[{1,1,n}, tRange, opts, NumberOfFrames->60
]/; NumberQ[N@n]

LimaconOfPascalGenerator2[
	dist_:1, tRange_List:{Automatic}, opts:OptionsPattern[]
] :=
	If[ tRange==={Automatic},
		ConchoidGenerator[ {Cos@#,Sin@#}&, {0,2 Pi}, {{-1,0}, dist}, opts],
		ConchoidGenerator[ {Cos@#,Sin@#}&, tRange, {{-1,0}, dist}, opts]
	] /; NumberQ[N@dist]

NephroidGenerator[ tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	EpitrochoidGenerator[{1,1/2,1/2}, tRange, opts, NumberOfFrames->60]

Options[ParabolaGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

ParabolaGenerator[tRange_List:{-3.2,3.2}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,tRangeList, pointList,
	circleGP, lineGP, pointGP, lastFrameGP
	},
	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin,tMax,tStep];
	pointList = {#,1/4 #^2}& /@ tRangeList;

	circleGP = { Hue[.18], Circle[#, 1 + Last@#]}&;
	lineGP = { Hue[.65], Line[{{0,1}, #, {First@#,-1}}] }&;
	pointGP =
		{	PointSize[.02], Hue[.75], Point[{0,1}],
			Hue[0], Point[#],
			Hue[.75], Point[{First@#,-1}]
		}&;

	lastFrameGP =
		{	Hue[.75], Line[{{2 tMin,-1}, {2 tMax,-1}}],
			Hue[0], Point/@pointList,
			circleGP@ Last@ pointList,
			lineGP@ Last@ pointList,
			pointGP@ Last@ pointList
		};

		Table[
				Graphics[
					{	Hue[.75], Line[{{2 tMin,-1}, {2 tMax,-1}}],
						Hue[0], Point /@ Take[pointList,t],
						circleGP@ pointList[[t]],
						lineGP@ pointList[[t]],
						pointGP@ pointList[[t]]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->1.2 {{tMin,tMax}, {-1, Max[1, 1/4 tMin^2, 1/4 tMax^2] }}
				],
			{t, Length@tRangeList}
		]
]

Options[QuadratrixOfHippiasGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

QuadratrixOfHippiasGenerator[tRange_List:{0.0001, Pi/2}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,
	tRangeList, trailPointList,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin, tMax, tStep];
	trailPointList = (2 #/Pi){Cot[#],1}& /@ tRangeList;

	staticGP = {
		Line[{{0,0},{1,0},{1,1},{0,1},{0,0}}],
		Hue[0,0,0], Circle[{0,0},1, {0, Pi/2}],
		PointSize[.02], Point[{0,0}]
	};

(* movingGP is a function for mapping into integers *)
	movingGP =
	{	Hue[.65], Thickness[.008],
		Line[{{0,trailPointList[[#,2]]},{1,trailPointList[[#,2]]}}],
		Hue[.4, 1, .5],
		Line[{ {0,0}, {Cos[#], Sin[#]}&@ tRangeList[[#]] } ],
		PointSize[.02], Point@ trailPointList[[#]] 
	}&;

	lastFrameGP = {staticGP,
		Hue[0], PointSize[.01], Point /@ trailPointList,
		movingGP@ Length@ tRangeList
	};

		Table[
				Graphics[
					{ staticGP,
					Hue[0], PointSize[.01],Point/@ Take[ trailPointList,t],
						movingGP[t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->{{-.1,1.1},{-.1,1.1}}
				],
			{t, Length@ tRangeList}
		]
]

RoseGenerator[ n_:3, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	HypotrochoidGenerator[
		{n/(n+1), (n-1)/(2 (n+1)), 1/2},
		tRange, opts, NumberOfFrames->61
	] /; (IntegerQ[n] && n>1)

RoseGenerator2[n_:3, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
	If[ tRange==={Automatic},
		ParaPolarGenerator[
			{Cos[n #], #}&,
			{0, If[OddQ[n], Pi, 2 Pi] },
			opts
		],
		ParaPolarGenerator[ {Cos[n #], #}&, tRange, opts]
	]/; (IntegerQ[n] && n>1)

Options[RoseGenerator2] =
	Join[
		{NumberOfFrames->61},
		Options[Graphics]
	];

RoseGenerator2[ n_:3, tRange_List:{Automatic}, opts:OptionsPattern[]] :=
Module[{nOfFrames,t,tMin,tMax,tStep,
	tRangeList, curveXValueFun, curveYValueFun,center,
	pointListGP,staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		Which[
			tRange === {Automatic},
				{ 0, If[OddQ@n, Pi, 2 Pi], N[(tMax-tMin)/(nOfFrames-1)]},
			Length@tRange === 3,
				N@tRange,
			True,
				N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1)}
		];

	tRangeList = Range[tMin, tMax, tStep];

	curveXValueFun = Compile[{t}, Cos[n t] Cos@t];
	curveYValueFun = Compile[{t}, Cos[n t] Sin@t];

	pointListGP =
		Point[{curveXValueFun@#, curveYValueFun@# }]& /@ tRangeList;

	staticGP = {Hue[.65,1,.7], Thickness[.007], Circle[{0,0}, 1]};
	movingGP =
		{	Hue[.65,1,.7], Thickness[.007],
			Line[{{0,0}, {Cos[tMin +tStep (#-1)], Sin[tMin +tStep(#-1)] }}],
			Hue[.3,.7,.8], Thickness[.01], Line[{{0,0}, pointListGP[[#,1]]}],
			Hue[.65,1,.7], PointSize[.04], Point[{0,0}],
			Hue[0], PointSize[.02], pointListGP[[#]]
		}&;

	lastFrameGP = {staticGP, movingGP@Length@tRangeList, Hue[0],pointListGP};

		Table[
				Graphics[
					{	staticGP, movingGP[t],
						Hue[0], Take[ pointListGP,t]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				Ticks->{{-1,-.5,.5,1},{-1,-.5,.5,1}},
				PlotRange->{{-1.05,1.05},{-1.05,1.05}}
				],
			{t, Length@tRangeList}
		]
] /; IntegerQ[n] && n>1

Options[SineGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

SineGenerator[ tRange_List:{0, 2 Pi}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,
	staticGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	staticGP = {
		Line[{{-1,-1},{-1,1}}],
		Line[{{0,-1},{0,1}}],
		Line[{{2 Pi,-1},{2 Pi,1}}],
		Line[{{-2,0},{2 Pi,0}}],
		Hue[0,0,0], Circle[{-1,0}, 1]
	};

(* movingGP is a function for mapping into reals (range of t) *)
	movingGP =
	{
		{	Hue[.35,1,.5], Thickness[.004],
			Line[{{Cos[#]-1, Sin[#]}, {0,Sin[#]} }]
		},
		{	Hue[.65,1,.8], Thickness[.004],
			Table[ Line[{{t,0},{t,-Sin[t-#]}}], {t,0, 2 Pi, 2 Pi/20}],
			Thickness[.008],
			First@ Plot[ -Sin[t-#], {t,0, 2 Pi}, PlotStyle->{ Hue[.65,1,.8]}]
		},
		{Hue[.65,1,.8], PointSize[.02], Point[{0,Sin[#]}] },
		{	Hue[0,1,.8], Thickness[.008],
			Line[{ {-1,0}, {Cos[#]-1, Sin[#]} }]
		},
		{Hue[0,1,.8], PointSize[.02], Point[{Cos[#]-1, Sin[#]}] }
	}&;

	lastFrameGP = {staticGP, movingGP@ tMax, PointSize[.02], Point[{-1,0}]};

		Table[
				Graphics[{staticGP, movingGP[t],
					PointSize[.02], Point[{-1,0}]
				},
				FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->False,
				PlotRange->{{-2,2 Pi},{-1.25,1.25}} 1.05],
			{t, tMin, tMax, tStep}
		]
]

Options[TractrixGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

TractrixGenerator[tRange_List:{0, 1.52}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,
	tRangeList,trailPointList,lineListGP, movingGP,lastFrameGP},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@ tRange===3,
			N@ tRange,
			N@ {First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList =
		Table[ t^(2/3),
			{t, tMin^(3/2), tMax^(3/2),
			(tMax^(3/2)-tMin^(3/2))/(tMax-tMin) tStep}
		];
	trailPointList = {Log[Sec@# + Tan@#] -Sin[#], Cos[#]}& /@ tRangeList;
	lineListGP =
		Line[{
			#, { First@# + Sqrt[1-Last[#]^2], 0}
		}]& /@ trailPointList;

	(*movingGP is a pure function on integers *)
	movingGP =
		{	Hue[0], Thickness[.009], lineListGP[[#]],
			Hue[.3,1,.5],Thickness[.011],
			Line[{ {0,0}, lineListGP[[#,1,2]] }],
			Hue[0], PointSize[.02], Point@ trailPointList[[#]],
			Hue[.3,1,.5], PointSize[.025], Point@ lineListGP[[#,1,2]]
		}&;
	lastFrameGP =
		{Hue[0,.5,1], lineListGP, Point/@ trailPointList,
			movingGP@ Length@ tRangeList};

		Table[
				Graphics[
					{Hue[0,.5,1], Take[ lineListGP,t],
						movingGP@ t
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				Ticks->{Range[10],{0,1,2}},
				PlotRange->{{-.1,trailPointList[[-1,1]] +1.1},{-.1,1.5}}
				],
			{t, Length@ tRangeList}
		]
]

Options[WitchOfAgnesiGenerator] =
	Join[
		{NumberOfFrames->30},
		Options[Graphics]
	];

WitchOfAgnesiGenerator[tRange_List:{-1.1,1.1}, opts:OptionsPattern[]] :=
Module[{nOfFrames,tMin,tMax,tStep,tRangeList, pointList,
	movingGP, staticGP, lastFrameGP
	},

	nOfFrames = OptionValue[ NumberOfFrames ];

	{tMin, tMax, tStep} =
		If[ Length@tRange===3,
			N@tRange,
			N@{First@tRange, Last@tRange, (tMax-tMin)/(nOfFrames-1) }
		];

	tRangeList = Range[tMin,tMax,tStep];
	pointList = {2 Tan@#, 2 Cos[#]^2}& /@ tRangeList;

	movingGP =
		{	Hue[.65], Line[{{0,0}, {2 Tan@#,2}, {2 Tan@#,0}}],
			Line[{{2 Tan@tMin-1, 2 Cos[#]^2}, {2 Tan@tMax+1, 2 Cos[#]^2}}],
			PointSize[.02], Point[{Sin[2 #], 1+Cos[2 #]}],
			Hue[0], Point[{2 Tan@#, 2 Cos[#]^2}]
		}&;

	staticGP =
		{	Hue[.75], Circle[{0,1},1],
			Line[{{2 Tan@tMin-1, 0},{ 2 Tan@tMax+1, 0}}],
			Line[{{2 Tan@tMin-1, 2},{ 2 Tan@tMax+1, 2}}],
			PointSize[.02], Point[{0,0}]
		};

	lastFrameGP =
			{	Hue[0], Point/@pointList,
				staticGP, movingGP[tMax]
			};

		Table[
				Graphics[
					{	Hue[0], Point /@ Take[pointList,t],
						staticGP, movingGP@tRangeList[[t]]
					}, FilterRules[ {opts}, Options[ Graphics ] ],
				AspectRatio->Automatic, Axes->True,
				PlotRange->{ {2 Tan@tMin-.1, 2 Tan@tMax+.1}, {-.5, 2.6 }}
				],
			{t, Length@tRangeList}
		]
]

End[]

EndPackage[]
