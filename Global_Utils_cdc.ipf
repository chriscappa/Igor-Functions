#pragma rtGlobals=1		// Use modern global access method.

// Chris's Global Utilities
// some are stolen from Scott Herndon
// Rules: 
//	1. Can be used by any other functions and can reference each other
//	2. Cannot use any functions that are not in the global utilities or built in to IGOR

// Functions
AllAxesFreePosZero(graphname)
AnnualCycle(timewave,datawavestr,[newwavestr,BW,maskwavestr]) // calculate an annual cycle
DiscriminatedAxisList( graphName, axis_Type )
DiscriminatedTraceList( graphName, axisName )
MeanWithNaN(NameStr,PrintVar) // will calculate the mean of a wave that includes nans; can just use wavestats instead
MakeNaNFromGraph() // will the point indicated by cursor A into a NaN; cannot be undone
MakeNaNSectionFromGraph() // will set all points in a trace on a graph between Csr(A) and Csr(B) to nan
PrimaryColors(graphName) // will set traces on a graph to a series of colors
PrimaryColors_RGB(graphName) // alternative to PrimaryColors
ScaleAllYAxes( graph_name) // will scale all y-axes to the range shown on the window
StackAllAxes(graphName, lo_twink, hi_twink) // to create stacked graphs
StackAllAxes2(graphName, lo_twink, hi_twink) // alternative to StackAllAxes
StackThreeAxes(graphName, space) // if you just have three axes to stack
HandyGraphButtons() // will put some handy buttons on the front graph
Julian2Seconds(TimeWave,Year,TimeName) // will convert a wave of julian day time to Igor time; does not deal well when years change
Julian2SecondsVar(TimeVar,Year) // will return a single time
Seconds2Julian_WrapAround(TimeWave) // Like Julian2Seconds, but deals with year changes better
Seconds2Julian(TimeWave) // converts Igor time to julian time
RunningMean(wvstr,npnts) // calculates a running mean around each point
RunningMean2(namewave,n) // an alternative to RunningMean
ExtractFromGraph(NameString,NewNameString)
MeanFromGraph() // prints the mean value between Csr(A) and Csr(B) from a wave on a graph
Diurnal3(timewave,datawavestr,[newwavestr,BW,maskwavestr]) // replaces Diurnal and Diurnal2; calculates diurnal profile
GraphDiel(thewavestr) // will graph a diurnal profile, with central value and +/- std deviation
Diurnal_RangeFromGraph(timewave,datawavestr,[newwavestr,BW]) // Will calculate diurnal profile for selected wave over range on graph
Diurnal_Range(timewave,datawavestr,[newwavestr,BW,start,stop]) // calculate diurnal profile for subset of entire wave
ExtractForBoxAndWhisker(refwavestr,datawavestr,start,stepsize,nsteps,[newwavestr]) // legacy
IsEven(var) // determines if a number is even
ConcatenateMatrix(w1, w2) // will concatenate two matrices, with the first wave changed
ConcatenateMatrixNewWave(w1, w2, wnew) // will concatenate two matrices, with a new wave made
SumLayers(matrix) // Will sum over layers in a 3D wave; use instead MatrixOp
WaveStatsWithMask(wv,mask,[quiet]) // Do wavestats, but with a data mask applied
BoxAndWhisker_Setup(refwavestr,datawavestr,start,stepsize,nsteps,[newwavestr,maskwavestr,minpoints,logscale]) // first step of creating box and whisker plot
BoxAndWhisker_AppendToGraph(wvstr,[LogScale,MyAxis]) // will append a box and whisker plot to an existing graph
BinAndFit(ywavestr,xwavestr,binwavestr,newwavestr,start,stepsize,nbins,[minpoints,logscale,zerointercept]) // fit data after binning
BinAndAveRatio(ywavestr,xwavestr,binwavestr,newwavestr,start,stepsize,nsteps,[minpoints,logscale,MaskWave]) // calculate binned averages of the ratio of two waves
GetMedian(mywave,[maskwave]) // calculate the median of a wave
NormalizeAlongLayers(wvstr,[newwvstr,LL]) // normalize a 3D wave along layers
StringSpace2Semicolon(text) // replace spaces or multiple semicolons with a single semicolon in a string
SetGraphPrefs() // sets CDC's prefered conditions for axes, size
KappaToGF(kappa,RH,DpDry) // convert kappa to a growth factor
GFtoKappa(GF,RH,DpDry) // convert a growth factor to a kappa value
WeightedAverage(w,ew) // calculate a weighted average
GraphWavesDifAxes // make a graph with multiple panels
AppendWavesDifAxes // append waves to a graph with multiple panels

/////////////////////////////////////////////////////////////////////////////////////
Menu "Macros"
	Submenu "Global Utils"
		"Handy Graph Buttons", HandyGraphButtons()
		"Set Mirror, Pos and StandOff",  MirrorStandPos()
		"MakeNanFromGraph /F9", /Q, MakeNaNFromGraph()
//		"GetStartNuc /F9", /Q,GetStartNuc()
//		"SetBadPoint /F9", /Q, setbad()
	End
END

/////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////

Function AllAxesFreePosZero(graphName)
	String graphName
	
	String this, ax_list = AxisLIst( graphName )
	Variable idex = 0, count = ItemsInList( ax_list )
	if( count > 0 )
		do	
			this = StringFromList( idex, ax_list )
			ModifyGraph freePos($this)=0
			idex += 1
		while( idex < count )
	endif
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////

Function/T DiscriminatedAxisList( graphName, axis_Type )
	string graphName, axis_Type
	String return_str =""
	Variable left = cmpstr( LowerStr(axis_type), "left" )
	Variable bottom = cmpstr( LowerStr(axis_type), "bottom" )
	Variable right = cmpstr( LowerStr(axis_type), "right" )
	Variable top = cmpstr( LowerStr(axis_type), "top" )
	if( (left!=0) %& (bottom!=0) %& (right!=0) %& (top!=0) )
		return return_str
	endif
	String list = AxisList( graphName )
	Variable idex = 0
	String info
	do
		info = StringByKey( "AXTYPE", AxisInfo( graphName, StringFromList( idex, list ) ))
		if( cmpstr( LowerStr(info), Lowerstr(axis_Type) ) != 0 )
			list = RemoveListItem( idex, list )
		else
			idex += 1
		endif
	while( idex < ItemsInList( list ))
	return_str = list
	return return_str
End

Function/T DiscriminatedTraceList( graphName, axisName )
	string graphName, axisName
	String return_str =""
	
	String list = TraceNameList( graphName, ";", 1 )
	Variable idex = 0
	String info
	do
		info = StringByKey( "YAXIS", TraceInfo( graphName, StringFromList( idex, list ), 0 ))
		if( cmpstr( LowerStr(info), Lowerstr(axisName) ) != 0 )
			list = RemoveListItem( idex, list )
		else
			idex += 1
		endif
	while( idex < ItemsInList( list ))
	return_str = list
	return return_str
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////

Function /T ListTracesOnAxis( axis_name, graph_name )
	String axis_name, graph_name
	
	String return_list = ""
	String list_of_traces = TraceNameList( graph_name, ";", 1 )
	Variable trace_index = 0, num_of_traces_to_check = ItemsInList( list_of_traces )
	String info_str, this_axis_name
	String name_of_trace
	do
		name_of_trace = StringFromList( trace_index, list_of_traces )
		info_str = TraceInfo( graph_name, name_of_trace , 0)
		this_axis_name = StringByKey( "YAXIS", info_str )
		if( cmpstr( this_axis_name, axis_name )== 0 )
			return_list = return_list + name_of_trace + ";"
		endif
		trace_index += 1
	while( trace_index < num_of_traces_to_check )
	
	return return_list
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MeanWithNaN(NameStr,PrintVar)
	String NameStr
	variable PrintVar
	WAVE NameW = $NameStr
	
	WaveStats/Q NameW
	variable NumPtsNoNaN = V_npnts
	variable NumPtsWithNaN = numpnts(NameW)
	
	variable m
	variable idex = 0
	variable Sum1 = 0
	FOR(m=1;m<=NumPtsWithNaN;m+=1)
		IF(numtype(NameW[m-1]) == 2)
			//Do Nothing
		ELSE
			Sum1 = Sum1 + NameW[m-1]
			idex = idex + 1
		ENDIF
	ENDFOR
	
	variable/G MeanNan
	MeanNan = Sum1/idex
	IF(printvar == 1)
		print MeanNan
	ENDIF
	
	return(MeanNan)
END FUNCTION

//////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION GeometricMean(thewave)
	wave thewave
	
	WaveStats/Q thewave
	variable NumPtsNoNaN = V_npnts
	variable NumPtsWithNaN = numpnts(thewave)
	
	variable m
	variable idex = 0
	variable Sum1 = 1
	FOR(m=1;m<=NumPtsWithNaN;m+=1)
		IF(numtype(thewave[m-1]) == 2)
			//Do Nothing
		ELSE
			Sum1 = Sum1 *thewave[m-1]
			idex = idex + 1
		ENDIF
	ENDFOR
	
	variable MeanNan
	MeanNan = Sum1^(1/idex)

	return(MeanNan)
END FUNCTION

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MakeNaNFromGraph()

	variable index = pcsr(A)
	WAVE wA = CsrWaveRef(A)
	wA[index] = NaN

END FUNCTION
//
FUNCTION MakeNaNSectionFromGraph()

	variable index1 = pcsr(A)
	variable index2 = pcsr(B)
	WAVE wA = CsrWaveRef(A)
	wA[index1,index2] = NaN

END FUNCTION

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function MirrorStandPos()
	Execute( "ModifyGraph mirror=2,standoff=0, freePos=0")
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function PrimaryColors(graphName,[special])
	String graphName
	variable special
	
	Wave/Z pc_red=pc_red
	Wave/Z pc_green=pc_green
	Wave/Z pc_blue=pc_blue
	Make/O/N=26 pc_red = {  65280, 52224, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
	Make/O/N=26 pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
	Make/O/N=26 pc_blue={    0, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}
	
	if(special == 1)
	
	pc_red = {512, 0, 48896, 65024, 65280}
	pc_green = {23296, 34816, 59904, 52992, 27136, 0}
	pc_blue = { 48384, 52224, 65280, 28672, 27136, 0}
	endif
	
	Variable num_axes, axis_dex = 0, dex
	Variable num_traces_this_axis, trace_dex = 0
	String axis_list, this_trace_list, this_trace, this_ax
	axis_list = DiscriminatedAxisList( graphName, "left" )
	num_axes = ItemsInList( axis_list )
	if( num_axes == 0 )
		return -1
	endif
	do
		this_trace_list = DiscriminatedTraceList( graphName, StringFromList( axis_dex, axis_list ) )
		num_traces_this_axis = ItemsInList( this_trace_list )
		if( num_traces_this_axis == 0 )
			// do nothing
		else
			trace_dex = 0
			do
				this_trace = StringFromList( trace_dex, this_trace_list )
				if( num_axes == 1 )
					dex = Floor( trace_dex/5 ) + Mod( trace_dex, 5 ) * 5
				else
					dex = axis_dex * 5 + trace_dex
				endif
				
				ModifyGraph rgb($this_trace)=(pc_red[dex], pc_green[dex], pc_blue[dex] )
				if( (trace_dex == 0) %&(num_axes != 1) )
					this_ax = StringFromList( axis_dex, axis_list )
					ModifyGraph axRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph tlblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph alblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph gridRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
				endif
				trace_dex += 1
			while( trace_dex < num_traces_this_axis )
		endif
		axis_dex += 1
	while( axis_dex < num_axes )
	
	KillWaves/Z pc_red, pc_green, pc_blue
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function ColorOptics(graphName,[special])
	String graphName
	variable special
	
	Wave/Z pc_red=pc_red
	Wave/Z pc_green=pc_green
	Wave/Z pc_blue=pc_blue
//	Make/O/N=26 pc_red = {  65280, 52224, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
//	Make/O/N=26 pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
//	Make/O/N=26 pc_blue={    0, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}
	Make/O/N=26 pc_red = { 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 65280, 52224, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
	Make/O/N=26 pc_green={   12800, 0, 0, 0, 0, 39168, 1, 2, 3, 4, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
	Make/O/N=26 pc_blue={   52224, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}

	if(special == 1)
	
	pc_red = {512, 0, 48896, 65024, 65280}
	pc_green = {23296, 34816, 59904, 52992, 27136, 0}
	pc_blue = { 48384, 52224, 65280, 28672, 27136, 0}
	endif
	
	Variable num_axes, axis_dex = 0, dex
	Variable num_traces_this_axis, trace_dex = 0
	String axis_list, this_trace_list, this_trace, this_ax
	axis_list = DiscriminatedAxisList( graphName, "left" )
	num_axes = ItemsInList( axis_list )
	if( num_axes == 0 )
		return -1
	endif
	dex = 0
	do
		this_trace_list = DiscriminatedTraceList( graphName, StringFromList( axis_dex, axis_list ) )
		num_traces_this_axis = ItemsInList( this_trace_list )

		if( num_traces_this_axis == 0 )
			// do nothing
		else
			trace_dex = 0
			do
				this_trace = StringFromList( trace_dex, this_trace_list )
				if( num_axes == 1 )
					dex = Floor( trace_dex/5 ) + Mod( trace_dex, 5 ) * 5
				else
					dex = axis_dex * 5 + trace_dex
				endif
				
				ModifyGraph rgb($this_trace)=(pc_red[dex], pc_green[dex], pc_blue[dex] )
				if( (trace_dex == 0) %&(num_axes != 1) )
					this_ax = StringFromList( axis_dex, axis_list )
					ModifyGraph axRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph tlblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph alblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph gridRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
				endif
				trace_dex += 1
			while( trace_dex < num_traces_this_axis )
		endif
		axis_dex += 1
	while( axis_dex < num_axes )
	
	KillWaves/Z pc_red, pc_green, pc_blue
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function PrimaryColors_RGB(graphName)
	String graphName
	
	Wave/Z pc_red=pc_red
	Wave/Z pc_green=pc_green
	Wave/Z pc_blue=pc_blue
//	Make/O/N=26 pc_red = {  65280, 52224, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
//	Make/O/N=26 pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
//	Make/O/N=26 pc_blue={    0, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}
	Make/O/N=26 pc_red = {  65280, 0, 0}
	Make/O/N=26 pc_green={   0, 50280, 0}
	Make/O/N=26 pc_blue={    0, 0, 60280}
	
	Variable num_axes, axis_dex = 0, dex
	Variable num_traces_this_axis, trace_dex = 0
	String axis_list, this_trace_list, this_trace, this_ax
	axis_list = DiscriminatedAxisList( graphName, "left" )
	num_axes = ItemsInList( axis_list )
	if( num_axes == 0 )
		return -1
	endif
	dex = 0
	do
		this_trace_list = DiscriminatedTraceList( graphName, StringFromList( axis_dex, axis_list ) )
		num_traces_this_axis = 1 //  ItemsInList( this_trace_list )
		if( num_traces_this_axis == 0 )
			// do nothing
		else
			trace_dex = 0
//			do
				this_trace = StringFromList( trace_dex, this_trace_list )
//				if( num_axes == 1 )
//					dex = Floor( trace_dex/5 ) + Mod( trace_dex, 5 ) * 5
//				else
//					dex = axis_dex * 5 + trace_dex
//				endif
				
				ModifyGraph rgb($this_trace)=(pc_red[dex], pc_green[dex], pc_blue[dex] )
				if( (trace_dex == 0) %&(num_axes != 1) )
					this_ax = StringFromList( axis_dex, axis_list )
					ModifyGraph axRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph tlblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph alblRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
					ModifyGraph gridRGB($this_ax)=(pc_red[dex], pc_green[dex], pc_blue[dex])
				endif
//				trace_dex += 1
//			while( trace_dex < num_traces_this_axis )
		endif
		axis_dex += 1
		dex += 1
	while( axis_dex < num_axes )
	
	KillWaves/Z pc_red, pc_green, pc_blue
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This next function will rescale all of the y-axes on a plot so that for a given
// x range, the 'autoscale' of over the present range

Function ScaleAllYAxes( graph_name)
	String graph_name
	
	Variable low_x, high_x
	Variable num_of_axes
	Variable num_of_traces
	
	Variable this_y_max, this_y_min // terminology, this refers to the current instance or trace
	Variable glob_y_max, glob_y_min // refers to the ultimate or universal max and min (all traces)
	String axes_list, this_axis_name, trace_list, this_trace_name
	Variable axis_index, trace_index
	
	String trace_x_info, trace_x_wave, trace_x_wave_df, x_axis_name, cmd
	Variable low_point, high_point
	
	// This is a bug, I think, do it later
	//GetAxis/W=$graph_name/Q bottom
	//low_x = V_min; high_x = V_max
	
	// part one determine how many left or right axes on plot
	axis_index = 0
	axes_list =DiscriminatedAxisList( graph_name, "left" ) + DiscriminatedAxisList( graph_name, "right" )
	
	num_of_axes = ItemsInList( axes_list )
	if( num_of_axes < 1 )
		//Print "Error in Global Utility: AutoYAllAxes -- too few axes"
		return -1
	endif
	//Wave wave_ref
	do
		// within part one step two cycle through all traces associated with that axis
		this_axis_name = StringFromList(axis_index, axes_list )
		trace_index = 0
		trace_list = ListTracesOnAxis( this_axis_name, graph_name )
		// determine y_u_max, and y_u_min
		GetAxis/W=$graph_name/Q $this_axis_name
		glob_y_min = V_max; glob_y_max = V_min // YES these are switched on purpose
	
		do
			this_trace_name = StringFromList( trace_index, trace_list )	
			trace_x_info = TraceInfo( graph_name, this_trace_name, 0 )
			trace_x_wave = StringByKey( "XWAVE", trace_x_info )
			x_axis_name = StringByKey( "XAXIS", trace_x_info )
			//printf "for %s :::%s on %s gives ", this_trace_name, x_axis_name, graph_name
			GetAxis/W=$graph_name/Q $x_axis_name
			low_x = V_min; high_x = V_max
			//printf "the min as %g and teh max as %g\r", low_x, high_x
			if( strlen( trace_x_wave )< 1 )
				// it is an x scaled wave and we get the information via
				WaveStats/Q /R=(low_x, high_x) TraceNameToWaveRef( graph_name, this_trace_name )
				
			else
				// it is an x - y pair and we need the infor mation from the x_wave 
				Wave wx = XwaveRefFromTrace( graph_name, this_trace_name )
				// with out of order pairs this bugs... I don't know how to fix it at the moment
				low_point = BinarySearch( wx, low_x )
				if( low_point == -1 )
					low_point = 0
				endif
				high_point = BinarySearch( wx, high_x )
				if( high_point == -2 )
					high_point = numpnts( wx )
				endif
				//printf "In Scale Y:  low_pnt:%g & low_x:%g t hi_pnt:%g and hi_x:%g", low_point, low_x, high_point, high_x
				WaveStats/Q /R=[low_point, high_point] TraceNameToWaveRef( graph_name, this_trace_name )
				//printf "   Result in %g, %g for Y's\r", V_min, V_max
			endif
			this_y_min = V_min; this_y_max = V_max
			if( this_y_min < glob_y_min )
				glob_y_min = this_y_min
			endif
			if( this_y_max > glob_y_max )
				glob_y_max = this_y_max
			endif
			trace_index += 1
		while( trace_index < ItemsInList( trace_list ) )
		//Print "Axis Name: "+this_axis_name+" y(max): "+num2str(glob_y_max)+" y(min): "+num2str(glob_y_min)
		
		// can we querry the axis name here
		// and if it is log scale, make sure we don't set to negative or zero
		String info_again = AxisInfo(graph_name, this_axis_name )
		String log_str = StringByKey( "log(x)", info_again, "=" )
		Variable log_val = str2num( log_str )
		if( log_val == 1 )
			if( glob_y_min < 0 )
				glob_y_min = 1e-9
			endif
		endif
		SetAxis/W=$graph_name/N=1 $this_axis_name glob_y_min, glob_y_max
		axis_index += 1
while( axis_index < num_of_axes  )
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Usage Notes for StackAllAxes
// StackAllAxes looks for all left axes and stacks them
// on graphName
// lo_twink and hi_twink are the number of hundredth to inwardly twink the middle axes by
Function StackAllAxes(graphName, lo_twink, hi_twink)
	String graphName
	Variable lo_twink, hi_twink
	
	String axes_list = AxisList( graphName ), this_ax, this_ax_type
	axes_list = RemoveFromList(  "top", axes_list)
	axes_list = RemoveFromList( "bottom", axes_list)
	Variable num_axes = ItemsInList( axes_list )
	
	if( num_axes == 0 )
		return -1
	endif
	Variable kdex = 0
	do
		this_ax = StringFromList( kdex, axes_list )
		this_ax_type = StringByKey( "AXTYPE", AxisInfo( graphName, this_ax ) )
		if( (cmpstr( lowerstr( this_ax_type), "top" ) == 0)      %|     (cmpstr( lowerstr( this_ax_type), "bottom" ) == 0) )
			axes_list = RemoveFromList( this_ax, axes_list )
		endif
		kdex += 1
	while( kdex < num_axes )
		 num_axes = ItemsInList( axes_list )
	Variable idex = 0, interval = 1/num_axes, low_bound, high_bound
	
	
	do
		this_ax = StringFromList( idex, axes_list )
		low_bound = idex * interval
		high_bound = (idex + 1) * interval
		if( (idex >= 0) %& (idex < num_axes + 1 ) )
	//		low_bound += lo_twink/100;	high_bound -= hi_twink/100 // changed by CD Cappa 07/30/10
			low_bound = low_bound; 		high_bound -= hi_twink/100
		endif
		if( low_bound < 0 )
			low_bound = 0
		endif
		if( high_bound > 1 )
			high_bound = 1 
		endif
		ModifyGraph axisEnab($this_ax)={low_bound,high_bound}
		idex += 1
	while( idex < num_axes )
	return 1
	
	End
	
//******************************************************************************
Function StackAllAxes2(graphName, lo_twink, hi_twink)
	String graphName
	Variable lo_twink, hi_twink
	
	String axes_list = AxisList( graphName ), this_ax, this_ax_type
	axes_list = RemoveFromList(  "top", axes_list)
	axes_list = RemoveFromList( "bottom", axes_list)
	Variable num_axes = ItemsInList( axes_list )
	
	if( num_axes == 0 )
		return -1
	endif
	Variable kdex = 0
	do
		this_ax = StringFromList( kdex, axes_list )
		this_ax_type = StringByKey( "AXTYPE", AxisInfo( graphName, this_ax ) )
		if( (cmpstr( lowerstr( this_ax_type), "top" ) == 0)      %|     (cmpstr( lowerstr( this_ax_type), "bottom" ) == 0) )
			axes_list = RemoveFromList( this_ax, axes_list )
		endif
		kdex += 1
	while( kdex < num_axes )
		 num_axes = ItemsInList( axes_list )
	Variable idex = 0, interval = 1/num_axes, low_bound, high_bound
	
	
	do
		this_ax = StringFromList( idex, axes_list )
		low_bound = idex * interval + 0.02
		high_bound = (idex + 1) * interval - 0.02
		if( (idex >= 0) %& (idex < num_axes + 1 ) )
	//		low_bound += lo_twink/100;	high_bound -= hi_twink/100 // changed by CD Cappa 07/30/10
			low_bound = low_bound; 		high_bound -= hi_twink/100
		endif
		if( idex == num_axes-1 )
			high_bound = 1
		endif
		if( idex == 0 )
			low_bound = 0
		endif
		if( low_bound < 0 )
			low_bound = 0
		endif
		if( high_bound > 1 )
			high_bound = 1 
		endif
		ModifyGraph axisEnab($this_ax)={low_bound,high_bound}
		idex += 1
	while( idex < num_axes )
	
	return 1
	
End

//******************************************************************************
Function StackThreeAxes(graphName, space)
	String graphName
	Variable space
	
	String axes_list = AxisList( graphName ), this_ax, this_ax_type
	axes_list = RemoveFromList(  "top", axes_list)
	axes_list = RemoveFromList( "bottom", axes_list)
	Variable num_axes = ItemsInList( axes_list )
	
	if( num_axes == 0 )
		return -1
	endif
	Variable kdex = 0
	do
		this_ax = StringFromList( kdex, axes_list )
		this_ax_type = StringByKey( "AXTYPE", AxisInfo( graphName, this_ax ) )
		if( (cmpstr( lowerstr( this_ax_type), "top" ) == 0)      %|     (cmpstr( lowerstr( this_ax_type), "bottom" ) == 0) )
			axes_list = RemoveFromList( this_ax, axes_list )
		endif
		kdex += 1
	while( kdex < num_axes )
		 num_axes = ItemsInList( axes_list )
	Variable idex = 0, interval = 1/num_axes, low_bound, high_bound
	
	variable counter
	do
		this_ax = StringFromList( idex, axes_list )
		low_bound = idex * interval + space/100
		high_bound = (idex + 1) * interval - space/100
		if(idex==0)
			low_bound = 0
		elseif(idex==num_axes-1)
			high_bound = 1
		endif
//		if( (idex >= 0) %& (idex < num_axes + 1 ) )
//	//		low_bound += lo_twink/100;	high_bound -= hi_twink/100 // changed by CD Cappa 07/30/10
//			low_bound = low_bound; 		high_bound -= hi_twink/100
//		endif
//		if( idex == num_axes-1 )
//			high_bound = 1
//		endif
//		if( idex == 0 )
//			low_bound = 0
//		endif
//		if( low_bound < 0 )
//			low_bound = 0
//		endif
//		if( high_bound > 1 )
//			high_bound = 1 
//		endif
		ModifyGraph axisEnab($this_ax)={low_bound,high_bound}
		ModifyGraph freePos($this_ax)={0,bottom},mirror($this_ax)=1,axThick($this_ax)=1.5,fSize($this_ax)=14
		idex += 1
	while( idex < num_axes )
		ModifyGraph axThick=1.5
		ModifyGraph fSize=14
	
	return 1
	
End
	
//*******************************************************	
	Function ColorTraceOrAxes( graphName, code)
	String graphName
	Variable code
	
	
	Make/O CCRed = 	{65280,        0,        0,        0, 65280, 52224, 30464, 0 }
	Make/O CCBlue = 	{       0, 52224, 15872, 65280, 43520,        0, 30464, 0 }
	Make/O CCGreen =	{       0,        0, 65280, 65280,        0, 20736, 30464, 0 }
	
	Make/O/N=16 red_chart, green_chart, blue_chart
	//					0		1		2		3		4		5			6		7		8		9		10		11		12		13		14		15
	red_chart = { 		0,	65280,	65280/4, 	0, 	65280/2, 	65280/2, 	0,	65280/2,		0,		0,	65280/2,	65280/2,		0,	65280/4,		0,		0}
	green_chart = { 		0,		0,	65280/2, 	0, 	65280/2, 		0, 	65280/2,		0,	65280/2,		0,	65280/2,		0,	65280/2,		0,		0,	65820/4}
	blue_chart = { 		0,		0,	65280/4, 65280, 	65280/2, 	65280, 	65280/2,		0,		0,	65280/2,		0,	65280/2,	65280/2,		0,	65280/4,		0}
	
	String axes_list = AxisList( graphName ), this_ax
	axes_list = RemoveFromList(  "top", axes_list)
	axes_list = RemoveFromList( "bottom", axes_list)
	
	String this_trace, trace_list = TraceNameList( graphName, ";", 1), trace_info, just_wave
	
	Variable instance, num_axes = ItemsInList( axes_list ), number_sign
	if( num_axes == 0 )
		return -1
	endif
	Variable idex = 0, jdex, correct_color = 0, axis_colored = 0
	do
		this_ax = StringFromList( idex, axes_list )
		jdex = 0
		do
			this_trace = StringFromList( jdex, trace_list )
			number_sign = strsearch( this_trace, "#", 0 )
			if( number_sign > 0 )
				just_wave = this_trace[0, number_sign - 1 ]
				instance = str2num( this_trace[number_sign+1, strlen( this_trace )] )
				trace_info = TraceInfo( graphName, just_wave, instance )
			else
				trace_info = TraceInfo( graphName, this_trace, 0 )
			endif
			if( cmpstr( StringbyKey( "YAXIS", trace_info ), this_ax ) == 0 )
				if( num_axes == 1 )
					correct_color = jdex
				else
					correct_color = idex
				endif
				ModifyGraph rgb($this_trace)=(CCRed[correct_color], CCBlue[correct_color], CCGreen[correct_color])
		
				if( (code %& 1) %& (axis_colored) )
					ModifyGraph axRGB($this_ax)=(CCRed[correct_color], CCBlue[correct_color], CCGreen[correct_color])
					ModifyGraph tlblRGB($this_ax)=(CCRed[correct_color], CCBlue[correct_color], CCGreen[correct_color])
					ModifyGraph alblRGB($this_ax)=(CCRed[correct_color], CCBlue[correct_color], CCGreen[correct_color])
					ModifyGraph gridRGB($this_ax)=(CCRed[correct_color], CCBlue[correct_color], CCGreen[correct_color])
					axis_colored = 1
				endif
			endif
			jdex += 1
		while( jdex < ItemsInList( trace_list ) )
		idex += 1 
	while( idex < num_axes )	
	
	KillWaves /Z CCRed, CCBlue, CCGreen, red_chart, blue_chart, green_chart

End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function UsefulButtons(ctrlName) : ButtonControl
	String ctrlName

Variable span, new_min, new_max

	String botaxlist = DiscriminatedAxisList( "", "bottom" )
	if( ItemsInList( botaxlist ) == 0 )
		print "In procedure, 'usefulButtons' - cannot determine a 'bottom' axis to use"
		return -1
	endif
	String anybotax = StringFromList( 0, botaxlist )
	if( cmpstr( ctrlName, "tog_left") == 0 )
		GetAxis /Q $anybotax
		span = V_max - V_min
		new_min = V_min - span/2
		new_max = V_max - span/2
		SetAxis $anybotax new_min, new_max
	endif
	if( cmpstr( ctrlName, "tog_right") == 0 )
		GetAxis /Q $anybotax
		span = V_max - V_min
		new_min = V_min + span/2
		new_max = V_max + span/2
		SetAxis $anybotax new_min, new_max
	endif
	if( cmpstr( ctrlName, "widen") == 0 )
		GetAxis /Q $anybotax
		span = V_max - V_min
		new_min = V_min - span/2
		new_max = V_max + span/2
		SetAxis $anybotax new_min, new_max
	endif
	if( cmpstr( ctrlName, "scale_y") == 0 )
		ScaleAllYAxes("")
	endif
End
//////////////////////////////////////////
Macro HandyGraphButtons()
	PauseUpdate; Silent 1		// building window...
	ShowInfo
	//ShowTools
	
	Variable width = 55, height = 21
	Variable line = 0, col_1 = 0 * width , col_2 = width + 5, col_3 = 2 * width + 5 , col_4 = 3 * width + 5
	
	ControlBar height + height * 1 // this 1 means there is only one line of buttons right now
	Button tog_left,pos={col_1, line},size={width, height},proc=UsefulButtons,title="<pan left<"
	Button tog_right,pos={col_2, line},size={width,height},proc=UsefulButtons,title=">pan right>"
	Button widen,pos={col_3, line},size={width,height},proc=UsefulButtons,title="<widen>"
	Button scale_y, pos={ col_4, line}, size={width, height}, proc=UsefulButtons, title = "Scale Y"
EndMacro

/////////////////////////////////////////////////////////////////////////////////////////////////////////

FUNCTION Julian2Seconds(TimeWave,Year,TimeName)
	WAVE TimeWave
	variable Year
	String TimeName
	
//	variable Yr,Mo,Da
	
	variable TimeMin = floor(TimeWave[0])
	variable SecsThisYear = (24*60*60)*(TimeMin-1)
	
	if(strlen(TimeName) == 0)
		TimeName = "TimeSecs"
	else
		// use name as entered
	endif
	duplicate/o/d TimeWave $TimeName
	wave TimeS = $TimeName
	TimeS = (TimeWave - TimeMin)*24*60*60//+ date2secs(Yr,Mo,Da-1)
	variable SecsSince1904 = date2secs(Year,1,1)
	TimeS = SecsSince1904+TimeS+SecsThisYear
	setscale/p d 0,1, "dat", TimeS
END FUNCTION

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

FUNCTION Julian2SecondsVar(TimeVar,Year)
	variable TimeVar
	variable Year
	
	variable Yr,Mo,Da
	variable/G TimeSecsVar
	
	variable Tmin = floor(TimeVar)
	variable SecsThisYear = (24*60*60)*(Tmin-1)
	TimeSecsVar = (TimeVar - Tmin)*24*60*60//+ date2secs(Yr,Mo,Da-1)
	variable SecsSince1904 = date2secs(Year,1,1)
	TimeSecsVar = SecsSince1904+TimeSecsVar+SecsThisYear
	Return TimeSecsVar
END FUNCTION

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Converts Time to Julian Day; Creates wave Time_JD
// If year changes, resets to zero
FUNCTION Seconds2Julian_WrapAround(TimeWave)
	WAVE TimeWave
	
	variable Yr,Mo,Da
	variable npnts = numpnts(TimeWave)
	variable m
	variable YearVar
	variable YearSecs
	variable SecsPerDay = 60*60*24
	
	make/o/T/n=(npnts) datestrwave
	make/o/d/n=(npnts) Time_JD = nan
		
	datestrwave = secs2date(TimeWave,-2)
	
	for(m=0;m<npnts;m+=1)
		sscanf datestrwave[m], "%i", yearvar
		YearSecs = date2secs(yearvar,1,1)
		Time_JD[m] = (TimeWave[m] - YearSecs)/SecsPerDay
	endfor

	KillWaves/Z datestrwave

END FUNCTION

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Converts Time to Julian Day; Creates wave Time_JD
// If year changes, keeps counting up
FUNCTION Seconds2Julian(TimeWave,[JDname])
	WAVE TimeWave
	string JDname
	
	if(paramisdefault(JDname))
		JDname = "Time_JD"
	endif
	
	variable npnts = numpnts(TimeWave)
	variable m
	variable YearVar
	variable YearSecs
	variable SecsPerDay = 60*60*24
	
	make/o/d/n=(npnts) $JDname = nan
	wave Time_JD = $JDname // Time_JD = nan
		
	string datestr = secs2date(TimeWave[0],-2)
	sscanf datestr, "%i", yearvar
	YearSecs = date2secs(yearvar,1,1)
	SecsPerDay = 60*60*24
	
	Time_JD = (TimeWave - YearSecs)/SecsPerDay

END FUNCTION

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION RunningMean(wvstr,npnts)
	string wvstr
	variable npnts	// number of points in average, must be odd
	variable npnts2 = (npnts-1)/2
	wave wv2avg = $wvstr
	variable nRows = numpnts(wv2avg)
	make/o/d/n=(nRows,4) $(wvstr+"_RM")
	wave newWv = $(wvstr+"_RM")

	variable avg
	variable i
	
	for(i=0;i<nRows;i+=1)
		wavestats/Q/R=[i-npnts2,i+npnts2] wv2avg
		newWv[i][0] = V_avg
		newWv[i][1] = V_avg+V_sdev
		newWv[i][2] = V_avg-V_sdev
		newWv[i][3] = V_npnts
	endfor
END FUNCTION

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION RunningMean_old(namewave,n,n2)
	string namewave
	variable n,n2
	// n = smoothing factor; n2 = offset
	WAVE NameW = $namewave
	string NewWaveStr = namewave + "_smth"
	Duplicate/O NameW $NewWaveStr
	WAVE NewWave = $NewWaveStr
	MAKE/O/N=(n)/D WaveForSmoothing

	variable V_npnts = numpnts($namewave)
	variable m, nn
	for(m=1;m<=V_npnts;m+=1)
		if((m-1) < n)
			// Do Nothing
		else
			WaveForSmoothing = NameW[p+(m-n-1)]
			MeanWithNaN("WaveForSmoothing",0)
			NVAR MeanNaN
			NewWave[m-1] = MeanNaN
			//WaveStats/Q/R=(m-n-1,m-1) NameW
			//NewWave[m-1] = V_avg
		endif
	endfor
	KillWaves WaveForSmoothing
END FUNCTION

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION RunningMean2(namewave,n)
	string namewave
	variable n
	// n = smoothing factor
	WAVE NameW = $namewave
	string NewWaveStr = namewave + "_smth"
	Duplicate/O NameW $NewWaveStr
	Duplicate/O NameW SmoothWave  // this is the wave that gets added to
	SmoothWave = 0
	WAVE NewWave = $NewWaveStr
	MAKE/O/N=(n)/D WaveForSmoothing

	variable V_npnts = n
	variable m, nn
	for(m=1;m<=V_npnts;m+=1)
		
		SmoothWave = SmoothWave + NameW[(x-n) + m]
		
	endfor
		NewWave = SmoothWave/n
		note/k NewWave "Smoothed using n = " + num2istr(n) + " from " + namewave
	KillWaves WaveForSmoothing, SmoothWave
END FUNCTION

///////////////////////////////////////////////////////////////////////////////////////////////////

FUNCTION ExtractFromGraph(NameString,NewNameString)
	String NameString, NewNameString
	
	
	IF(strlen(NameString) == 0)  // Determine whether to use wave from graph cursors or as entered
		NameString = csrwave(A)
	ELSE
		// Use NameString as entered
	ENDIF
	
	IF(WaveExists($NameString) == 1)  // Determine whether wave exists before continuing
		// continue
	ELSE
		print "Wave does not exist, cannot execute procedure"
		abort
	ENDIF
	
	IF(strlen(NewNameString) == 0) // Determine what to append to wavename
		NewNameString = "Ex"
	ELSE
		// Use NewNameString as entered
	ENDIF
	
	WAVE WaveToExtract = $NameString
	NewNameString = NameString + "_" + NewNameString
	//WAVE ExtractedWave = $NewNameString
	
	string stringA = csrwave(A)
	string stringB = csrwave(B)
	variable stringCompare = cmpstr(stringA,stringB)
	IF(stringCompare == 0) // Determine whether cursors are on the same wave in the graph
		// continue
	ELSE
		print "WARNING: Cursors are not on the same wave"
		abort
	ENDIF
	
	Duplicate/O/R=[pcsr(A),pcsr(B)] WaveToExtract $NewNameString
	print "Created wave " + NewNameString
	
END

///////////////////////////////////////////////////////////////////////////////////
Function MeanFromGraph()

//Need to make sure that cursors are on same trace!
WAVE wA = CsrWaveRef(A)
WAVE wB = CsrWaveRef(B)
WAVE xwave = CsrWaveRef(A)
wavestats/q/r=(xcsr(A),xcsr(B)) xwave
string result = "Mean is  " +num2str(V_avg)+ "  +/- "+num2str(V_sdev)
print result
return V_avg
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION GraphMatrix(Wave2GraphY, Wave2GraphX, ZDim,NewOrAppend)
	wave Wave2GraphY, Wave2GraphX
	variable ZDim, NewOrAppend
	
	variable NumColumns = dimsize(Wave2GraphY,1)
	variable m

IF(NewOrAppend == 1)	
	for(m=1;m<=NumColumns;m+=1)
		if(m==1)
			appendtograph Wave2GraphY[][m-1][ZDim] vs Wave2GraphX
		else
			appendtograph Wave2GraphY[][m-1][ZDim] vs Wave2GraphX
		endif
	endfor
ELSE
	for(m=1;m<=NumColumns;m+=1)
		if(m==1)
			Display Wave2GraphY[][m-1][ZDim] vs Wave2GraphX
		else
			appendtograph Wave2GraphY[][m-1][ZDim] vs Wave2GraphX
		endif
	endfor
ENDIF
	
END FUNCTION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Requires time in Julian Day for proper operation...see "Seconds2Julian"
function Diurnal(timewave,datawave,newwavestr)
	wave timewave, datawave
	string newwavestr
	
	variable nhours = 24
	variable m
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days
		seconds2julian(timewave)
		wave TimeJD = Time_JD
		kill = 1
	else
		wave TimeJD = timewave
	endif
	
	make/o/d/n=(nhours) $newwavestr
	wave newwave = $newwavestr
	string str = newwavestr + "_sem"
	make/o/d/n=(nhours) $str
	wave newwave_sem = $str
	
	make/o/d/n=(nhours) Time_diurnal
	
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
	extract datawave, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
	newwave[m] = meanwithnan("diurnalwave",0)
	wavestats/q diurnalwave
	newwave[m] = V_avg
	newwave_sem[m] = V_sem
	Time_diurnal[m] = (Time_i+time_f)*24/2
	endfor
	
	killwaves/z diurnalwave
	if(kill==1)
		killwaves/z TimeJD
	endif
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Time wave can be seconds or julian day
// can create matrix for box and whisker plotting; replicates are stored in columns
function Diurnal2(timewave,datawavestr,[newwavestr,BW])
	wave timewave		// time wave (seconds or julian day)
	string datawavestr	// data wave
	string newwavestr		// new output wave name; optional
	variable BW			// 1 = create wave for box and whisker plotting
	
	string BWwavestr
	
	if(ParamIsDefault(newwavestr)==1)
		newwavestr = datawavestr + "_diel"
		BWwavestr = datawavestr + "_dielBW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	if(ParamIsDefault(BW)==1)
		BW = 0
	endif
	
	wave datawave = $datawavestr
	
	if(waveexists(datawave)==0 || waveexists(timewave)==0)
		abort "one of your waves don't exist."
	endif
	
	variable nhours = 24
	variable m
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days; make julian day wave and kill at end
		seconds2julian(timewave)
		wave TimeJD = Time_JD
		kill = 1
	else
		wave TimeJD = timewave
	endif
	
	make/o/d/n=(nhours,3) $newwavestr
	wave newwave = $newwavestr
	setscale/P x, 0.5, (24/nhours), newwave
	
	variable np
	variable nt = numpnts(timewave)
	variable npmax=0
	variable nd
		
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		wavestats/q diurnalwave
		newwave[m][0] = V_avg
		newwave[m][1] = V_sem
		newwave[m][2] = V_sdev
		np = numpnts(diurnalwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		np = numpnts(diurnalwave)
		if(BW==1)
			if(m==0)
				make/o/d/n=(npmax,nhours) $BWwavestr = nan
				wave DielResultForBW = $BWwavestr
			endif
			DielResultForBW[0,np-1][m] = diurnalwave[p]
		endif		
	endfor
		
	killwaves/z diurnalwave, timepts_temp
	if(kill==1)
		killwaves/z TimeJD
	endif
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Time wave can be seconds or julian day
// can create matrix for box and whisker plotting; replicates are stored in columns
function Diurnal3(timewave,datawavestr,[newwavestr,BW,maskwavestr])
	wave timewave		// time wave (seconds or julian day)
	string datawavestr	// data wave
	string newwavestr		// new output wave name; optional
	variable BW			// 1 = create wave for box and whisker plotting
	string maskwavestr	// 1 = good data
	
	string BWwavestr
	
	if(ParamIsDefault(newwavestr)==1)
		newwavestr = datawavestr + "_diel"
		BWwavestr = datawavestr + "_dielBW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	if(ParamIsDefault(BW)==1)
		BW = 0
	endif
	if(ParamIsDefault(maskwavestr)==1)
		make/o/d/n=(numpnts(timewave))/FREE masktemp = 1
		wave maskwave = masktemp
	else
		wave maskwave = $maskwavestr
	endif
	
	wave datawave = $datawavestr
	
	if(waveexists(datawave)==0 || waveexists(timewave)==0)
		abort "one of your waves don't exist."
	endif
	
	variable nhours = 24
	variable m
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days; make julian day wave and kill at end
		seconds2julian(timewave)
		wave TimeJD = Time_JD
		kill = 1
	else
		wave TimeJD = timewave
	endif
	
	make/o/d/n=(nhours,3) $newwavestr
	wave newwave = $newwavestr
	setscale/P x, 0.5, (24/nhours), newwave
	
	duplicate/o/d datawave datawavetemp
	datawavetemp = maskwave==1 ? datawavetemp : nan
	
	variable np
	variable nt = numpnts(timewave)
	variable npmax=0
	variable nd
		
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawavetemp, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		wavestats/q diurnalwave
		newwave[m][0] = V_avg
		newwave[m][1] = V_avg+V_sdev//V_sem
		newwave[m][2] = V_avg-V_sdev//V_sdev
		np = numpnts(diurnalwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawavetemp, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		np = numpnts(diurnalwave)
		if(BW==1)
			if(m==0)
				make/o/d/n=(npmax,nhours) $BWwavestr = nan
				wave DielResultForBW = $BWwavestr
			endif
			DielResultForBW[0,np-1][m] = diurnalwave[p]
		endif		
	endfor
		
	if(ParamIsDefault(maskwavestr)==0)
		note/k newwave "Data have been masked by " + maskwavestr
	endif

	
	killwaves/z diurnalwave, timepts_temp
	if(kill==1)
		killwaves/z TimeJD
	endif
End

//************************************************************************************
Function GraphDiel(thewavestr)
	string thewavestr
	
	wave thewave = $thewavestr
	
	display thewave[][1], thewave[][2], thewave[][0]
	ModifyGraph mode($thewavestr)=7,hbFill($thewavestr)=2,toMode($thewavestr)=1
	ModifyGraph rgb($(thewavestr+"#2"))=(0,0,0),lsize($(thewavestr+"#2"))=2
	setaxis left 0,*
	setaxis bottom 0,24
	ModifyGraph mirror=1,fSize=14,axThick=1.5,standoff=0
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Time wave can be seconds or julian day
// can create matrix for box and whisker plotting; replicates are stored in columns
function Diurnal_RangeFromGraph(timewave,datawavestr,[newwavestr,BW])
	wave timewave		// time wave (seconds or julian day)
	string datawavestr	// data wave
	string newwavestr		// new output wave name; optional
	variable BW			// 1 = create wave for box and whisker plotting
	
	string BWwavestr
	
	if(ParamIsDefault(newwavestr)==1)
		newwavestr = datawavestr + "_diel"
		BWwavestr = datawavestr + "_dielBW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	if(ParamIsDefault(BW)==1)
		BW = 0
	endif
	
	wave datawave = $datawavestr
	// select range from graph
	make/o/d/n=(numpnts(datawave))/FREE datawave_temp
	datawave_temp = datawave
	datawave_temp[,min(pcsr(A),pcsr(B))] = nan
	datawave_temp[max(pcsr(A),pcsr(B)),] = nan
	
	if(waveexists(datawave)==0 || waveexists(timewave)==0)
		abort "one of your waves don't exist."
	endif
	
	variable nhours = 24
	variable m
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days; make julian day wave and kill at end
		seconds2julian(timewave)
		wave TimeJD = Time_JD
		kill = 1
	else
		wave TimeJD = timewave
	endif
	
	make/o/d/n=(nhours,3) $newwavestr
	wave newwave = $newwavestr
	setscale/P x, 0.5, (24/nhours), newwave
	
	variable np
	variable nt = numpnts(timewave)
	variable npmax=0
	variable nd
		
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave_temp, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		wavestats/q diurnalwave
		newwave[m][0] = V_avg
		newwave[m][1] = V_sem
		newwave[m][2] = V_sdev
		np = numpnts(diurnalwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		np = numpnts(diurnalwave)
		if(BW==1)
			if(m==0)
				make/o/d/n=(npmax,nhours) $BWwavestr = nan
				wave DielResultForBW = $BWwavestr
			endif
			DielResultForBW[0,np-1][m] = diurnalwave[p]
		endif		
	endfor
	
	string range_info = "Selected between points " + num2str(min(pcsr(A),pcsr(B))) + " and " + num2str(max(pcsr(A),pcsr(B)))
	note/k newwave range_info	
	
	killwaves/z diurnalwave, timepts_temp
	if(kill==1)
		killwaves/z TimeJD
	endif
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Time wave can be seconds or julian day
// can create matrix for box and whisker plotting; replicates are stored in columns
function Diurnal_Range(timewave,datawavestr,[newwavestr,BW,start,stop])
	wave timewave		// time wave (seconds or julian day)
	string datawavestr	// data wave
	string newwavestr		// new output wave name; optional
	variable BW			// 1 = create wave for box and whisker plotting
	variable start			// index to start
	variable stop			// index to stop
	
	string BWwavestr
	
	if(ParamIsDefault(newwavestr)==1)
		newwavestr = datawavestr + "_diel"
		BWwavestr = datawavestr + "_dielBW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	if(ParamIsDefault(BW)==1)
		BW = 0
	endif
	
	wave datawave = $datawavestr
	
	if(waveexists(datawave)==0 || waveexists(timewave)==0)
		abort "one of your waves don't exist."
	endif
	
	make/o/d/n=(numpnts(datawave))/FREE datawave_temp
	// select range from graph
	if(ParamIsDefault(start)==1)
		start = 0
	endif
	if(ParamIsDefault(stop)==1)
		stop = numpnts(datawave)
	endif
	datawave_temp = datawave
	datawave_temp[,start] = nan
	datawave_temp[stop,] = nan
	
	variable nhours = 24
	variable m
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days; make julian day wave and kill at end
		seconds2julian(timewave)
		wave TimeJD = Time_JD
		kill = 1
	else
		wave TimeJD = timewave
	endif
	
	make/o/d/n=(nhours,3) $newwavestr
	wave newwave = $newwavestr
	setscale/P x, 0.5, (24/nhours), newwave
	
	variable np
	variable nt = numpnts(timewave)
	variable npmax=0
	variable nd
		
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave_temp, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		wavestats/q diurnalwave
		newwave[m][0] = V_avg
		newwave[m][1] = V_sem
		newwave[m][2] = V_sdev
		np = numpnts(diurnalwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	for(m=0;m<nhours;m+=1)
		time_i = m/24
		time_f = (m+1)/24
		extract datawave, diurnalwave, (TimeJD[x] - floor(TimeJD[x])) >= time_i && (TimeJD[x] - floor(TimeJD[x])) < time_f
		np = numpnts(diurnalwave)
		if(BW==1)
			if(m==0)
				make/o/d/n=(npmax,nhours) $BWwavestr = nan
				wave DielResultForBW = $BWwavestr
			endif
			DielResultForBW[0,np-1][m] = diurnalwave[p]
		endif		
	endfor
	
	string range_info = "Selected between points " + num2str(start) + " and " + num2str(stop)
	note/k newwave range_info	
	
	killwaves/z diurnalwave, timepts_temp
	if(kill==1)
		killwaves/z TimeJD
	endif
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate a matrix that can be used for making Box and Whisker plots
// Replicates are stored as COLUMNS
function ExtractForBoxAndWhisker(refwavestr,datawavestr,start,stepsize,nsteps,[newwavestr])
	string refwavestr			// name of the wave that contains the values to which you want to bin your data
	string datawavestr		// data wave
	variable start				// start of binning range. 
	variable stepsize			// bin width
	variable nsteps			// number of bins
	string newwavestr			// new name for output wave; optional
	
	string BWwavestr
	string BWxstr
	
	if(ParamIsDefault(newwavestr)==1)
		BWwavestr = datawavestr + "_BW"
		BWxstr = datawavestr + "_BWxwave"
	else
		BWwavestr = newwavestr + "_BW"
		BWxstr = newwavestr + "_BWx"
	endif
	
	wave refwave = $refwavestr
	wave datawave = $datawavestr
	
	variable m
	variable np
	variable npmax=0
	variable low, high
		
	for(m=0;m<nsteps;m+=1)
		low = start + m*stepsize
		high = start + (m+1)*stepsize
		extract datawave, tempwave, refwave >= low && RefWave < high
		np = numpnts(tempwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	make/o/d/n=(npmax,nsteps) $BWwavestr = nan
	wave BWresult = $BWwavestr
	make/o/d/n=(nsteps) $BWxstr
	wave BWx = $BWxstr

	for(m=0;m<nsteps;m+=1)
		low = start + m*stepsize
		high = start + (m+1)*stepsize
		extract datawave, tempwave, refwave >= low && RefWave < high
		np = numpnts(tempwave)
		BWresult[0,np-1][m] = tempwave[p]
		BWx[m] = low + 0.5*stepsize
	endfor
		
	killwaves/z tempwave
End

///////////////////////////////////////////////////////////////////////
Function IsEven(var)
	variable var
	
	variable var2
	var2 = var/2
	if(floor(var2) == var2)
		return 1
	else
		return 0
	endif
End

///////////////////////////////////////////////////////////////////////
Function IsDivisible(var,divisor)
	variable var
	variable divisor
	
	variable var2
	var2 = var/divisor
	if(floor(var2) == var2)
		return 1
	else
		return 0
	endif
End

//*******************************************************************
Function ConcatenateMatrix(w1, w2)
//		Tacks the contents of w2 on the to end of w1.
//		If w1 does not exist, it is created.
//		This is designed for 2D waves.
	String w1, w2
	
	Variable rows1, rows2
	Variable cols1, cols2

	if (Exists(w1) == 0)
		Duplicate $w2, $w1
	else
		String wInfo=WaveInfo($w2, 0)
		Variable WType=NumberByKey("NUMTYPE", wInfo)
		rows1 = dimsize($w1,0)//numpnts($w1)
		rows2 = dimsize($w2,0)//numpnts($w2)
		cols1 = dimsize($w1,1)
		cols2 = dimsize($w2,1)
		if(cols1!=cols2)
			abort "2D waves must have same number of columns."
		endif
		Redimension/N=((rows1 + rows2),-1) $w1  // preserve number of columns
		if (WType)				// Numeric wave
			Wave/C/D ww1=$w1
			Wave/C/D ww2=$w2
			ww1[rows1, ][] = ww2[p-rows1][q]
		else						// Text wave
			Wave/T tw1=$w1
			Wave/T tw2=$w2
			tw1[rows1, ][] = tw2[p-rows1][q]
		endif
	endif
End

//*******************************************************************
Function ConcatenateMatrixNewWave(w1, w2, wnew)
//		Tacks the contents of w2 on the to end of w1.
//		If w1 does not exist, it is created.
//		This is designed for 2D waves.
	String w1, w2, wnew
	
	Variable rows1, rows2
	Variable cols1, cols2

	if (Exists(w1) == 0)
		Duplicate $w2, $w1
	else
		String wInfo=WaveInfo($w1, 0)
		Variable WType=NumberByKey("NUMTYPE", wInfo)
		rows1 = dimsize($w1,0)//numpnts($w1)
		rows2 = dimsize($w2,0)//numpnts($w2)
		cols1 = dimsize($w1,1)
		cols2 = dimsize($w2,1)
		if(cols1!=cols2)
			abort "2D waves must have same number of columns."
		endif
		if (WType)				// Numeric wave
			if(Wtype==2)
				make/o/n=((rows1+rows2),cols1) $wnew = nan
			elseif(Wtype==1)
				make/o/c/n=((rows1+rows2),cols1) $wnew = nan
			else//(Wtype==4)
				make/o/d/n=((rows1+rows2),cols1) $wnew = nan
			endif
			Wave/C/D wwnew=$wnew
			Wave/C/D ww1=$w1
			Wave/C/D ww2=$w2
			wwnew[,rows1-1][] = ww1[p][q]
			wwnew[rows1, ][] = ww2[p-rows1][q]
		else						// Text wave
			make/o/t/n=((rows1+rows2),cols1) $wnew = ""
			wave/T twnew=$wnew
			Wave/T tw1=$w1
			Wave/T tw2=$w2
			twnew[,rows1-1][] = tw1[p][q]
			twnew[rows1, ][] = tw2[p-rows1][q]
		endif
	endif
End

//**************************************************************
Function SumLayers(matrix)
	// will sum a 3D matrix by the layers index
	// result is a 2D matrix
	wave matrix
	
	variable nRows = dimsize(matrix,0)
	variable nCols = dimsize(matrix,1)
	variable nLayers = dimsize(matrix,2)
	
	make/o/d/n=(nRows,nCols) SummedMatrix
	variable t1, t2, t3
	variable i
	for(i=0;i<nLayers;i+=1)
		SummedMatrix += matrix[p][q][i]
	endfor
	
//	return(SummedMatrix)
End

//********************************************************************************
Function WaveStatsWithMask(wv,mask,[quiet])
	wave wv		// the wave to get properties of
	wave mask	// your mask wave, where 1 = good data
	variable quiet	// 0 (default) = print results; 1 = quiet
	
	if(ParamIsDefault(quiet))
		quiet = 0
	endif
	if(numpnts(wv)!=numpnts(mask))
		abort "Your waves have to have the same number of points. Try again."
	endif
	if(waveexists(wv)==0)
		abort "Your wave doesn't exist."
	elseif(waveexists(mask)==0)
		abort "Your mask doesn't exist."
	endif
	make/o/d/n=(numpnts(wv))/FREE WaveWithMask
	WaveWithMask = mask == 1 ? wv : nan
	
	if(quiet==0)
		wavestats WaveWithMask
	else
		wavestats/q WaveWithMask
	endif	
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate a matrix that can be used for making Box and Whisker plots
function BoxAndWhisker_Setup(refwavestr,datawavestr,start,stepsize,nsteps,[newwavestr,maskwavestr,minpoints,logscale])
	string refwavestr			// name of the wave that contains the values to which you want to bin your data
	string datawavestr		// data wave
	variable start				// start of binning range. 
	variable stepsize			// bin width
	variable nsteps			// number of bins
	string newwavestr			// new name for output wave; optional
	string maskwavestr		// name of wave to use for data masking; optional
	variable MinPoints		// minimum number of points that must exist in a given bin range
	variable LogScale			// treat stepsize according to log scaling
	
	if(ParamIsDefault(MinPoints))
		MinPoints = 5
	endif
	if(ParamIsDefault(LogScale))
		LogScale = 0
	endif
	
	string BWwavestr
	string BWxstr
	string BWdf
	string cdf = getdatafolder(2)
	
	if(ParamIsDefault(newwavestr)==1)
		BWwavestr = datawavestr + "_BW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	
	wave refwave = $refwavestr
	wave datawave = $datawavestr
	variable nRows = numpnts(datawave)

	if(waveexists(refwave)==0 || waveexists(datawave)==0)
		abort "your waves don't exist. choose wisely."
	endif

	if(ParamIsDefault(newwavestr)==1)
		BWdf = "BW_"+datawavestr
		if(datafolderexists(BWdf)==0)
			newdatafolder $(BWdf)
		endif
	else
		BWdf = "BW_"+newwavestr
		if(datafolderexists(BWdf)==0)
			newdatafolder $(BWdf)
		endif
	endif
	
	if(ParamIsDefault(maskwavestr))
		make/o/d/n=(nRows)/FREE mask4BW = 1
		wave mask = mask4BW
		maskwavestr = ""
	else
		wave mask = $maskwavestr
	endif	
	
	make/o/d/n=(nRows)/FREE datawavetemp = nan, refwavetemp=nan
	datawavetemp = mask==1 ? datawave : nan
	refwavetemp = mask==1 ? refwave : nan 
	
	variable m
	variable np
	variable npmax=0
	variable low, high
	variable median, meanval, firstq_upper, secondq_upper, firstq_lower, secondq_lower, xleft, xright
		
	setdatafolder $BWdf
	make/o/d/n=(nsteps,13) $BWwavestr = nan
	wave BWresult = $BWwavestr
	setscale/P x, start, stepsize, BWresult
	if(stringmatch(maskwavestr,"")==0)
		note/k BWresult "Binned wave " + datawavestr + " according to " + refwavestr + " with maskwave = " + maskwavestr
	else
		note/k BWresult "Binned wave " + datawavestr + " according to " + refwavestr + " with no masking"
	endif
	
	note BWresult "Bin parameters are: " + num2str(start) + "; " + num2str(stepsize) + "; " + num2str(nsteps) + "; " + num2str(logscale)
	
	make/o/d/n=(nSteps*3,3) BWfill = nan
	variable stepsize_init = stepsize
	make/o/d/n=(nsteps) RangeMidpt

//	if(LogScale==1)
//		make/o/d/n=(nsteps) RangeStart, RangeStop
////		make/o/d/n=(nsteps) RangeMidpt
//		for(m=0;m<nsteps;m+=1)
//			RangeStart[m] = 10^(log(start + m*stepsize - stepsize/2)) // start + stepsize*10^(x*LogScale)
//			RangeStop[m] = 10^(log(start + m*stepsize + stepsize/2)) // start + stepsize*10^((x+1)*LogScale) // start + 10^(stepsize+(x+1)*LogScale)//(x+1)*10^(stepsize+(x+1)*LogScale)
//		endfor
//		RangeMidpt = (RangeStop+RangeStart)/2
	if(logscale!=0)

		make/o/d/n=(nsteps) RangeStart, RangeStop
		// old method
//		RangeStart = 10^(log(start + (x)^logscale*stepsize - stepsize/2))
//		RangeStop = RangeStart[x+1]
//		RangeStop[numpnts(RangeStop)-1] = 10^(log(start + (numpnts(RangeStop))^logscale*stepsize - stepsize/2))
//		RangeMidpt = (RangeStop+RangeStart)/2
		wavestats/q refwavetemp
		variable minx = V_min
		variable maxx = V_max
		variable step = ((log(V_max)-log(V_min))/(nsteps-1)
		RangeStart = 10^(log(V_min) + step*x)
		RangeStop = 10^(log(V_min) + step*(x+1)*0.9999)
		RangeMidpt = (RangeStop+RangeStart)/2
	endif


	for(m=0;m<nsteps;m+=1)
		if(LogScale==0)
			low = start + m*stepsize - stepsize/2
			high = start + (m+1)*stepsize - stepsize/2
			RangeMidpt[m] = (low+high)/2
		else
			low = RangeStart[m]
			high = RangeStop[m]
			stepsize = RangeStop[m]-RangeStart[m]
		endif
		extract datawavetemp, tempwave, refwave >= low && RefWave < high && numtype(datawavetemp)!=2
		if(numpnts(tempwave)>=MinPoints)
			sort tempwave, tempwave	// sort smallest to largest
			wavestats/q tempwave 
			redimension/N=(V_npnts) tempwave
			
			BWresult[m][0] = tempwave[V_npnts/2]		// median
			BWresult[m][1] = V_avg					// mean
			BWresult[m][2] = tempwave[V_npnts/4]		// lower quartile
			BWresult[m][3] = BWresult[m][0] - BWresult[m][2]		// lower quartile error bars
			BWresult[m][4] = tempwave[3*V_npnts/4]	// upper quartile
			BWresult[m][5] = tempwave[3*V_npnts/4] - BWresult[m][0]	// upper quartile error bars
			BWresult[m][6] = BWresult[m][2] - tempwave[V_npnts/10]	// 10th percentile error bars
			BWresult[m][7] = tempwave[9*V_npnts/10] - BWresult[m][4]	// 90th percentile error bars
			BWresult[m][8] = 0.8*(high-low)/2		// lower x value
			BWresult[m][9] = 0.8*(high-low)/2		// cannot recall why this is the same as 8
			BWresult[m][10] = V_npnts			// number of points
			BWresult[m][11] = V_sdev				// standard deviation
			BWresult[m][12] = V_sem				// standard deviation of the mean
			BWfill[m*3,m*3+1][0] = BWresult[m][4]
			BWfill[m*3,m*3+1][1] = BWresult[m][2]
			BWfill[m*3][2] = low+stepsize/2-BWresult[m][8]*0.95
			BWfill[m*3+1][2] = low+stepsize/2+BWresult[m][8]*0.95
		else
			BWresult[m][] = nan
			BWfill[m*3,m*3+1][] = nan
		endif
	endfor

	killwaves/z tempwave
	setdatafolder $cdf
End

//***************************************************************************************
Function BoxAndWhisker_AppendToGraph(wvstr,[LogScale,MyAxis,NumOnGraph,appendMeanOnly,appendMedOnly])
	string wvstr
	variable LogScale
	string MyAxis
	variable NumOnGraph // how many B&W graphs are there on your graph so far?
	variable appendMeanOnly // 1 = append mean values
	variable appendMedOnly // 1 = append median values only
	
	if(ParamIsDefault(LogScale))
		LogScale = 0
	endif
	variable DifferentAxis
	if(ParamIsDefault(MyAxis))
		DifferentAxis = 0
		MyAxis = "left"
	else
		DifferentAxis = 1
	endif
	if(ParamIsDefault(NumOnGraph))
		NumOnGraph=0
	endif
	
	string fldrBW = "BW_"+wvstr
	String fldrSav0= GetDataFolder(1)
	String checkit = stringfromlist(itemsinlist(fldrSav0,":")-1,fldrSav0,":")
	string BWstr
	if(stringmatch(checkit,fldrBW)!=1)
		// you are in the wrong folder go down one level
		BWstr = wvstr+"_BW"
		wvstr = ":"+fldrBW+":"+wvstr+"_BW"
		wave BWwave = $(wvstr)//":"+fldrBW+":"+wvstr+"_BW")
		wave BWfill = $(":"+fldrBW+":BWfill")
		wave RangeMidpt = $(":"+fldrBW+":RangeMidpt")
	else
		// you are in the right folder
		BWstr = wvstr+"_BW"
		wvstr  = wvstr+"_BW"
		wave BWwave = $(wvstr)//+"_BW")
		wave BWfill
		wave RangeMidpt = $(":"+fldrBW+":RangeMidpt")
	endif
	if(waveexists(BWwave)==0)
		abort "your wave don't exist, fool."
	endif
	
	NumOnGraph = NumOnGraph*5
	string BWstr1 = BWstr+"#1"
	string BWstr2 = BWstr+"#2"
	string BWstr3 = BWstr+"#3"
	string BWstr4 = BWstr+"#4"
	string BWfillstr = "BWfill"
	string BWfillstr1 = "BWfill#1"
	
	if(appendMeanOnly!=1 && appendMedOnly!=1)
		// append full box and whisker
		if(DifferentAxis==0)
			AppendToGraph BWfill[][0], BWfill[][1] vs BWfill[][2]
		else
			AppendToGraph/l=$MyAxis BWfill[][0], BWfill[][1] vs BWfill[][2]
		endif
		ModifyGraph mode($BWfillstr)=7,usePlusRGB($BWfillstr)=1,hbFill($BWfillstr)=2,toMode($BWfillstr)=1
		ModifyGraph rgb($BWfillstr)=(0,0,0),plusRGB($BWfillstr)=(65535,65535,65535)
		ModifyGraph rgb($BWfillstr1)=(0,0,0)
		
		if(DifferentAxis==0)
			if(LogScale==0)
				AppendToGraph BWwave[*][0],BWwave[*][0],BWwave[*][2],BWwave[*][4],BWwave[*][1]
			else
				wave RangeMidpt = $(":"+fldrBW+":RangeMidpt")
				AppendToGraph BWwave[*][0],BWwave[*][0],BWwave[*][2],BWwave[*][4],BWwave[*][1] vs RangeMidpt
			endif
		else
			if(LogScale==0)
				AppendToGraph/R=$MyAxis BWwave[*][0],BWwave[*][0],BWwave[*][2],BWwave[*][4],BWwave[*][1]
			else
				wave RangeMidpt = $(":"+fldrBW+":RangeMidpt")
				AppendToGraph/R=$MyAxis BWwave[*][0],BWwave[*][0],BWwave[*][2],BWwave[*][4],BWwave[*][1] vs RangeMidpt
			endif
		endif
		ModifyGraph mode($BWstr)=2,mode($BWstr1)=2,mode($BWstr2)=2
		ModifyGraph mode($BWstr3)=2,mode($BWstr4)=3
		ModifyGraph marker($BWstr)=16,marker($BWstr4)=16
		ModifyGraph lSize($BWstr)=0
		ModifyGraph rgb($BWstr)=(0,0,0),rgb($BWstr1)=(0,0,0)
		ModifyGraph rgb($BWstr2)=(0,0,0),rgb($BWstr3)=(0,0,0),rgb($BWstr4)=(0,0,0)
		ErrorBars/T=2/L=0/Y=2 $BWstr BOX,wave=(BWwave[*][8],BWwave[*][9]),wave=(BWwave[*][5],BWwave[*][3])
		ErrorBars/T=0/L=1.5 $BWstr1 X,wave=(BWwave[*][8],BWwave[*][9])
		ErrorBars/T=1.5/L=1.5 $BWstr2 Y,wave=(,BWwave[*][6])
		ErrorBars/T=1.5/L=1.5 $BWstr3 Y,wave=(BWwave[*][7],)
	elseif(appendMeanOnly==1)
//		if(differentaxis==0)
//			AppendToGraph BWwave[*][1] vs RangeMidpt
//		else
			AppendToGraph/L=$MyAxis BWwave[*][1] vs RangeMidpt
//		endif
		if(numongraph>0)
			ModifyGraph mode($Bwstr + "#"+num2istr(numongraph))=0
		else
			ModifyGraph mode($BWstr)=0
		endif
	elseif(appendMedOnly==1)
		AppendToGraph/L=$myaxis BWwave[*][0] vs RangeMidpt
	else
		//oops
	endif
	
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate linear fits (slopes and intercepts) from two waves, but where the fit range
// is determined by a third wave to which the data are "binned"
function BinAndFit(ywavestr,xwavestr,binwavestr,newwavestr,start,stepsize,nbins,[minpoints,logscale,zerointercept])
	string ywavestr			// y data wave
	string xwavestr			// x data wave
	string binwavestr			// reference wave to use for binning
	string newwavestr			// new name for output wave; optional
	variable start				// start of binning range. 
	variable stepsize			// bin width
	variable nbins			// number of bins
	variable MinPoints		// minimum number of points that must exist in a given bin range
	variable LogScale			// treat stepsize according to log scaling
	variable zeroIntercept		// 1 = fix intercept to zero; 0 = free fit
	
	if(ParamIsDefault(MinPoints))
		MinPoints = 5
	endif
	if(ParamIsDefault(LogScale))
		LogScale = 0
	endif
	if(ParamIsDefault(zeroIntercept))
		zerointercept = 0
	endif
	
	string mywavestr
	string BWxstr
	string BWdf
	string cdf = getdatafolder(2)
	
	mywavestr = newwavestr + "_slope"
	
	wave binwave = $binwavestr
	wave ywave = $ywavestr
	wave xwave = $xwavestr
	variable nRows = numpnts(ywave)

	if(waveexists(ywave)==0 || waveexists(xwave)==0 || waveexists(binwave)==0)
		abort "your waves don't exist. choose wisely."
	endif
	
	// create mask wave (which will be updated as things progress)
	make/o/d/n=(nRows)/FREE mask4fit = 1
	wave mask = mask4fit
	make/o/d/n=(nRows,nbins) mask_byBin
	make/o/d/n=(nRows) yerr, xerr
	
	
	yerr = abs(0.2*ywave)
//	yerr = yerr < 1 ? 1 : yerr
	xerr = abs(0.2*xwave)
//	xerr = xerr < 0.2 ? 0.2 : xerr
		
	variable m
	variable np
	variable npmax=0
	variable low, high
	variable median, meanval, firstq_upper, secondq_upper, firstq_lower, secondq_lower, xleft, xright
		
	make/o/d/n=(nbins,5) $mywavestr = nan
	wave mywave = $mywavestr
	setscale/P x, start+stepsize/2, stepsize, mywave	// set the scale
	note mywave "Fit of " + ywavestr + " vs. " + xwavestr + " and binned according to " + binwavestr
	
	variable stepsize_init = stepsize
	
	if(LogScale!=0)
		make/o/d/n=(nbins)/FREE RangeStart, RangeStop
		make/o/d/n=(nbins)/FREE RangeMidpt
		RangeStart = start + x*10^(stepsize+x*LogScale)
		RangeStop = start + (x+1)*10^(stepsize+(x+1)*LogScale)
		RangeMidpt = (RangeStop+RangeStart)/2
	endif

	for(m=0;m<nbins;m+=1)
		if(LogScale==0)
			low = start + m*stepsize - stepsize/2
			high = start + (m+1)*stepsize - stepsize/2
		else
			low = RangeStart[m]
			high = RangeStop[m]
			stepsize = RangeStop[m]-RangeStart[m]
		endif
		
		mask4fit = (binwave>=low && binwave<high) ? 1 : 0
		mask_bybin[][m] = ywave[p]*mask4fit[p]
		wavestats/q mask4fit
		if(V_sum >= minpoints)
			if(zerointercept==0)
				CurveFit/Q/ODR=2 line  ywave /X=xwave /M=mask4fit /D /I=1 ///W=yerr /XW=xerr
			else
				K1 = 1
				CurveFit/H="01"/NTHR=0/Q/ODR=2 line  ywave /X=xwave /M=mask4fit /D /I=1 /W=yerr /XW=xerr
				note mywave, "Intercept fixed at zero"
			endif
			wave W_coef, W_sigma
			mywave[m][0] = W_coef[1] // slope
			mywave[m][1] = W_coef[0] // intercept
			mywave[m][2] = W_sigma[1] // fit error in slope
			mywave[m][3] = V_sum // number of points in fit
			mywave[m][4] = V_r2 // R^2
		else
			mywave[m][] = nan
		endif
	endfor
	
	mask_bybin = mask_bybin==0 ? nan : mask_bybin

	wave fitwave = $("fit_"+ywavestr)
	killwaves/z tempwave, W_coef, W_sigma, fitwave
	setdatafolder $cdf
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate binned averages of the ratio of two waves, but where the averages 
// are calculated over narrow ranges determined by a third wave to which the data are "binned"
function BinAndAveRatio(ywavestr,xwavestr,binwavestr,newwavestr,start,stepsize,nsteps,[minpoints,logscale,MaskWave])
	string ywavestr			// y data wave
	string xwavestr			// x data wave
	string binwavestr			// reference wave to use for binning
	string newwavestr			// new name for output wave; optional
	variable start				// start of binning range. 
	variable stepsize			// bin width
	variable nsteps			// number of bins
	variable MinPoints		// minimum number of points that must exist in a given bin range
	variable LogScale			// treat stepsize according to log scaling
	wave MaskWave		// optional mask wave (1=use data, 0 = do not use)
	
	if(ParamIsDefault(MinPoints))
		MinPoints = 5
	endif
	if(ParamIsDefault(LogScale))
		LogScale = 0
	endif
	variable useMask
	if(ParamIsDefault(MaskWave))
		useMask = 0
	else
		useMask = 1
	endif
		
	string mywavestr
	string BWxstr
	string BWdf
	string cdf = getdatafolder(2)
	
	mywavestr = newwavestr + "_ratio"
	
	wave binwave = $binwavestr
	wave ywave = $ywavestr
	wave xwave = $xwavestr
	variable nRows = numpnts(ywave)

	if(waveexists(ywave)==0 || waveexists(xwave)==0 || waveexists(binwave)==0)
		abort "your waves don't exist. choose wisely."
	endif
	
	// create mask wave (which will be updated as things progress)
	make/o/d/n=(nRows)/FREE mask4fit = 1
	wave mask = mask4fit
	make/o/d/n=(nRows) tempwave
		
	variable m
	variable np
	variable npmax=0
	variable low, high
	variable median, meanval, firstq_upper, secondq_upper, firstq_lower, secondq_lower, xleft, xright
		
	make/o/d/n=(nsteps,5) $mywavestr = nan
	wave mywave = $mywavestr
	setscale/P x, start+stepsize/2, stepsize, mywave	// set the scale
	note mywave "Fit of " + ywavestr + " vs. " + xwavestr + " and binned according to " + binwavestr
	
	variable stepsize_init = stepsize
	
	if(LogScale!=0)
		make/o/d/n=(nsteps)/FREE RangeStart, RangeStop
		make/o/d/n=(nsteps)/FREE RangeMidpt
		RangeStart = start + x*10^(stepsize+x*LogScale)
		RangeStop = start + (x+1)*10^(stepsize+(x+1)*LogScale)
		RangeMidpt = (RangeStop+RangeStart)/2
	endif

	variable npnts
	
	for(m=0;m<nsteps;m+=1)
		if(LogScale==0)
			low = start + m*stepsize - stepsize/2
			high = start + (m+1)*stepsize - stepsize/2
		else
			low = RangeStart[m]
			high = RangeStop[m]
			stepsize = RangeStop[m]-RangeStart[m]
		endif
		
		mask4fit = (binwave>=low && binwave<high) ? 1 : 0
		if(useMask==1)
			mask4fit*=maskwave
			note mywave "Mask applied"
		endif
		wavestats/q mask4fit
		npnts = V_sum
		if(V_sum >= minpoints)
			tempwave = (mask4fit==1) ? (ywave/xwave) : nan
			wavestats/q tempwave
			mywave[m][0] = V_avg // average
			mywave[m][1] = V_sdev // standard deviation
			mywave[m][2] = V_sem // standard error of the mean
			mywave[m][3] = npnts // number of points in fit
			mywave[m][4] = getmedian(tempwave) // median
		else
			mywave[m][] = nan
		endif
	endfor

	wave fitwave = $("fit_"+ywavestr)
	killwaves/z tempwave, W_coef, W_sigma, fitwave
	setdatafolder $cdf
End
//*************************************************************************************
Function GetMedian(mywave,[maskwave])
	// function to determine the median from an input wave
	// wave my contain nans
	// can also use a "masking" wave (0 = exclude, 1 = include)
	wave mywave
	wave maskwave
	
	variable doMask
	if(ParamIsDefault(maskwave))
		doMask = 0
	else
		doMask = 1
		if(numpnts(mywave)!=numpnts(maskwave))
			abort "length of your data wave and mask wave do not match."
		endif
	endif
	
	if(doMask==0)
		extract mywave, tempmedwave, numtype(mywave)!=2 // must have points
	else
		extract mywave, tempmedwave, numtype(mywave)!=2 && maskwave == 1
	endif
	sort tempmedwave, tempmedwave // sort smallest to largest
	variable MedianValue
	wavestats/q tempmedwave
	MedianValue = tempmedwave[round(V_npnts/2)]
	
	killwaves/z tempmedwave
	
	return MedianValue

End

//********************************************************************************
// Function to take a 3D wave and normalize each beam (i.e. row/column combination through layers)
// Normalizes to max value
Function NormalizeAlongLayers(wvstr,[newwvstr,LL])
	string wvstr	// wave to use
	string newwvstr	// enter here if you want to use a non-standard name
					// wave name is wvstr_norm by default (where wvstr is the name of the wave used)
	variable LL // lower limit to consider data as good
	
	if(ParamIsDefault(LL))
		LL = 0
	endif
	
	if(ParamIsDefault(NewWvStr))
		newwvstr = wvstr+"_norm"
	endif
	
	wave wv = $wvstr
	if(waveexists(wv)==0)
		abort "no such wave exists you fool."
	endif
	
	variable nRows = dimsize(wv,0)
	variable nCols = dimsize(wv,1)
	variable nLayers = dimsize(wv,2)
	make/o/d/n=(nRows,nCols,nLayers) $(newwvstr)
	wave normwv = $(newwvstr)
	note normwv "Normalized along layers from wave " + wvstr
	make/o/d/n=(nLayers)/FREE current
	variable i, j
	
	for(i=0;i<nRows;i+=1)
		for(j=0;j<nCols;j+=1)
			current = wv[i][j][p]
			wavestats/q current
			if(V_max>LL)
				normwv[i][j][] = current[r]/V_max
			else
				normwv[i][j][] = 0
			endif
		endfor
	endfor
end

//*************************************************************************
Function/T StringSpace2Semicolon(text)
	// take a string and remove any spaces, replace them with semi-colons, remove any double
	// or triple semicolons, and remove if there is a semi-colon as the first element. 
	string text
	
	variable length
	text = replacestring(" ",text,";")
	text = replacestring(";;;;;",text,";")
	text = replacestring(";;;;",text,";")
	text = replacestring(";;;",text,";")
	text = replacestring(";;",text,";")
	if(stringmatch(text[0],";"))
		length = strlen(text)
		text = text[1,length]
	endif
	
	return text
End

//*********************************************************************************
Function SetGraphPrefs()
	ModifyGraph tick=2,mirror=1,fSize=14,axThick=1.5,standoff=0
End

//******************************************************************************
Function KappaToGF(kappa,RH,DpDry)
// function to convert kappa values to growth factors
// based on Eqn. 11 in Petters and Kreidenweis (2007)
	variable kappa
	variable RH	// in percent
	variable DpDry	// dry diameter, nm
	
	variable surfacetension = 0.072 // J m^-2
	variable MWh2o = 0.018 // kg/mol
	variable IdealGas = 8.314 // J/(mol.K)
	variable Temperature = 298 // K
	variable densityH2O = 1000 // kg/m^3 
	variable Aparam = 1e9*(4*surfacetension*MWh2o)/(IdealGas*Temperature*densityH2O) // nm
	
	variable maxiters = 1000
	variable step = 0.0005
	variable GF = 1
	
	variable left, right

	if(numtype(Kappa)==2)
		GF = nan
	else	
		do
			left = 1e-2*RH/exp(Aparam/(DpDry*GF))
			right = (GF^3-1)/(GF^3-(1-kappa))
			
			if(GF>3)
				break
			elseif(right>left)
				break
			endif
			GF += step
		while(1)
	endif
	
	Return GF
End

//***********************************************************************
Function GFtoKappa(GF,RH,DpDry)
// function to convert GF values to kappa values
// based on Eqn. 11 in Petters and Kreidenweis (2007)
	variable GF
	variable RH	// in percent
	variable DpDry	// dry diameter, nm
	
	variable surfacetension = 0.072 // J m^-2
	variable MWh2o = 0.018 // kg/mol
	variable IdealGas = 8.314 // J/(mol.K)
	variable Temperature = 298 // K
	variable densityH2O = 1000 // kg/m^3 
	variable Aparam = 1e9*(4*surfacetension*MWh2o)/(IdealGas*Temperature*densityH2O) // nm
	
	variable B = exp(Aparam/(DpDry*GF))
	variable Kappa
	
	Kappa = (B/(RH*1e-2))*(GF^3-1)-GF^3+1
	
	Return Kappa
End

Function ppb_to_ugm3(ppb,MW)
	variable ppb
	variable MW
	
	variable conc
	conc = 1000000*ppb*MW*2.5E+25*0.000000001/6.022E+23
	return(conc)
End

//***************************************************************************
Function FindOutliers(thewave, s_thresh)
	// in progress. end effects not good yet. 
	wave thewave
	variable s_thresh // 1, 2 or 3 sigma
	
	variable npnts = numpnts(thewave)
	variable i
	variable npnts_per_group = 10 // decide how many to include
	
	make/o/d/n=(npnts) my_stdwave = nan, my_runningmean = nan, my_mask = 0
	make/o/d/n=(npnts_per_group+1) my_tempwave
	
	for(i=0;i<npnts;i+=1)
	
		if(i<npnts_per_group/2)
			my_tempwave = thewave[x+i]
		elseif(i > npnts-npnts_per_group)
			my_tempwave = thewave[x+i-npnts_per_group]
		else
			my_tempwave = thewave[x+i-npnts_per_group/2]
		endif
		
		my_tempwave[npnts_per_group/2] = nan
		wavestats/q my_tempwave
		my_stdwave[i] = v_sdev
		my_runningmean[i] = v_avg
		if(thewave[i] < my_runningmean[i] - s_thresh*my_stdwave[i])
			my_mask[i] = 1
		elseif(thewave[i] > my_runningmean[i] + s_thresh*my_stdwave[i])
			my_mask[i] = 1
		endif
	
	endfor
	
	
End


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate diurnal profiles of a time-series
// Will ignore NaN values
// Time wave can be seconds or julian day
// can create matrix for box and whisker plotting; replicates are stored in columns
function AnnualCycle(timewave,datawavestr,[newwavestr,BW,maskwavestr])
	wave timewave		// time wave (seconds only)
	string datawavestr	// data wave
	string newwavestr		// new output wave name; optional
	variable BW			// 1 = create wave for box and whisker plotting
	string maskwavestr	// 1 = good data
	
	string BWwavestr
	
	if(ParamIsDefault(newwavestr)==1)
		newwavestr = datawavestr + "_ann"
		BWwavestr = datawavestr + "_annBW"
	else
		BWwavestr = newwavestr + "_BW"
	endif
	if(ParamIsDefault(BW)==1)
		BW = 0
	endif
	if(ParamIsDefault(maskwavestr)==1)
		make/o/d/n=(numpnts(timewave))/FREE masktemp = 1
		wave maskwave = masktemp
	else
		wave maskwave = $maskwavestr
	endif
	
	wave datawave = $datawavestr
	
	if(waveexists(datawave)==0 || waveexists(timewave)==0)
		abort "one of your waves don't exist."
	endif
	
	variable nmonths = 12
	variable m, i
	variable time_i, time_f
	variable kill = 0
	
	if(timewave[0] > 366)
		// time is not in julian days; make julian day wave and kill at end
//		seconds2julian(timewave)
//		wave TimeJD = Time_JD
//		kill = 1
	else
		abort "time is in julian day, must be in igor time" // wave TimeJD = timewave
	endif
	
	make/o/d/n=(nmonths,3) $newwavestr
	wave newwave = $newwavestr
	setscale/P x, 1, (nmonths/nmonths), newwave
	
	duplicate/o/d datawave datawavetemp
	datawavetemp = maskwave==1 ? datawavetemp : nan
	
	variable np
	variable nt = numpnts(timewave)
	variable npmax=0
	variable nd
	
	make/o/t/n=(numpnts(timewave)) DateWaveStr
	make/o/d/n=(numpnts(timewave)) month_wv, year_wv
	string mydatestr
	
	DateWaveStr = secs2date(timewave,0)
	for(i=0;i<numpnts(timewave);i+=1)
		mydatestr = DateWaveStr[i]
		month_wv[i] = str2num(stringfromlist(0,mydatestr,"/"))
		year_wv[i] = str2num(stringfromlist(2,mydatestr,"/"))
	endfor
	
	variable current_month
	for(m=0;m<nmonths;m+=1)
		current_month = m + 1
		extract datawavetemp, diurnalwave, month_wv == current_month
		wavestats/q diurnalwave
		newwave[m][0] = V_avg
		newwave[m][1] = V_avg+V_sdev//V_sem
		newwave[m][2] = V_avg-V_sdev//V_sdev
		np = numpnts(diurnalwave)
		if(np>npmax)
			npmax = np
		endif
	endfor
	
	variable nyears
	wavestats/q year_wv
	variable current_year = V_min
	nYears = V_max - V_min + 1
	make/o/d/n=(nMonths,nYears,2) $BWwavestr = nan
	wave BWwave = $BWwavestr
	setscale/P x, 1, 1, BWwave
	setscale/P y, current_year, 1, BWwave

	
	for(i=0;i<nYears;i+=1)
		for(m=0;m<nmonths;m+=1)
			current_month = m + 1
			extract datawavetemp, diurnalwave, month_wv == current_month && year_wv == current_year
			if(numpnts(diurnalwave)>0)
				wavestats/q diurnalwave
				BWwave[m][i][0] = V_avg
				BWwave[m][i][1] = V_sem
			endif
		endfor
		current_year += 1
	endfor
		
	if(ParamIsDefault(maskwavestr)==0)
		note/k newwave "Data have been masked by " + maskwavestr
	endif

	
	killwaves/z diurnalwave, timepts_temp
//	if(kill==1)
//		killwaves/z TimeJD
//	endif
End

//************************************************************************************
Function GraphAnnualTrends(thewavestr)
	// will graph annual trends for the BW matrix generated from AnnualCycle()
	string thewavestr
	
	wave thewave = $thewavestr
	variable nMonths = dimsize(thewave,0)
	variable i
	
	for(i=0;i<nMonths;i+=1)
		if(i==0)
			display thewave[i][][0]
		else
			appendtograph thewave[i][][0]
		endif
	endfor
	setaxis left 0,*
	ModifyGraph mirror=1,fSize=14,axThick=1.5,standoff=0
end

//*************************************************************
Function InterpolateWaveCutEnds(OrigWave,ResultWaveStr)
// Uses NOAA interpolate wave function
// linear interpolation between points
// cuts off ends, so no extrapolation
	wave OrigWave
	string ResultWaveStr
	
	variable npnts = numpnts(origwave)
	Make/o/d/n=(npnts) $ResultWaveStr
	wave ResultWave = $ResultWaveStr
	note/k ResultWave "Interpolated from" 
	
	wavestats/q OrigWave
	if(V_npnts==0)
		abort
	endif
	
	variable i
	variable start, stop
	i = 0
	do
		if(numtype(origwave[i])!=2)
			start = i
			break
		else
			i += 1
		endif
	while(1)
	
	i = npnts-1
	do
		if(numtype(origwave[i])!=2)
			stop = i
			break
		else
			i -= 1
		endif
	while(1)
				
	InterpolateWave(ResultWave,OrigWave)
	resultwave[0,start-1] = nan
	resultwave[stop+1,] = nan
End

//*********************************************************************
Function Make_interp_between_pts(waveA,x0,x1)
	wave waveA
	
	variable x0, x1
	variable y0, y1, m
//	string waveAStr=CsrWave(A)
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "the cursors need to be on the same wave - Aborting from Make_interp_between_cursors"
	endif
	
//	WAVE waveA = CsrWaveRef(A)

	x0=pcsr(A)
	x1=pcsr(B)
	
	y0=waveA[x0]
	y1=waveA[x1]

	m = (y1-y0)/(x1-x0)

	if ( (numtype(y0)!=0) || (numtype(y1)!=0) )
		Abort "The cursors need to be on non-nan values - Aborting from Make_interp_between_cursors"
	endif
	
	if (x0 == x1 )
		Abort "The cursors need to be on different values - Aborting from Make_interp_between_cursors"
	endif
	if (x0>x1)
		abort "Cursor B needs to be ahead of cursor A."
	endif

	waveA[x0, x1]= y0 + m*(p-x0)

//	Print "// You inserted interpolated values in points ",x0," through ",x1," in wave "+waveAStr

End

//**************************************************************
Function/c WeightedAverage(w,ew)
	// calculate the weighted average of a wave
     wave w // wave with values to average
     wave ew // wave with errors
 
     variable av, su
     variable s1, s2
 
     duplicate/FREE w wtmp, ewtmp
     ewtmp = 1/ew^2
     wtmp = ewtmp*w
     wavestats/q wtmp
     s1 = V_sum
     wavestats/q ewtmp
     s2 = V_sum
     av = s1/s2 
     su = sqrt(numpnts(w))/sqrt(s2) 
 
     variable/c avw = cmplx(av,su) 
    return avw
 
end

//******************************************************************************
function GraphWavesDifAxes(ywvstr,xwvstr,ncols)
// make a multi-panel graph, where each wave goes on a different panel
	string ywvstr // semi-colon delimited string of waves to graph
	string xwvstr //  semi-colon delimited string of waves to graph
	variable ncols // number of columns
	
	variable nwvs = itemsinlist(ywvstr,";") // number of waves
	variable nrows = ceil(nwvs/ncols) // number of rows
	
	variable col_offset = 0.05
	variable col_len = (1-col_offset)*1/ncols
	variable col_space = ((1-col_len*ncols)/3)
	variable col_len2 = 1/ncols
	
	col_len = (1-col_offset)/ncols
	col_space = (1-col_len*ncols)/3
	col_space = 0.06
	col_len2 = col_len-col_space
	
	variable row_len
	if(nrows>1)
		row_len = (1-0.1/(nrows/2))*1/nrows
	else
		row_len = 1
	endif
	variable row_space = (1-row_len*nrows)/3
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	string xstr, ystr
	variable i, j
	variable wvcount = 0
	xstr = stringfromlist(0,xwvstr,";")
	xstr = "stuff"
	string graphname = stringfromlist(0,ywvstr,";") + "_v_" + xstr// stringfromlist(0,xwvstr,";")
	display/N=$graphname as graphname
	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			wvcount +=1
			if(wvcount > nwvs)
				abort
			endif
			yaxstr = "y" + num2istr(j+1 + ncols*i)
			xaxstr = "x" + num2istr(j+1 + ncols*i)
			xstr = stringfromlist(i*ncols+j,xwvstr,";")
			ystr = stringfromlist(i*ncols+j,ywvstr,";")
			wave xwv = $xstr
			wave ywv = $ystr
			appendtograph/l=$yaxstr/b=$xaxstr ywv vs xwv
//			print 1-row_len*(i+1)-row_space*i,1-(row_len)*i
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1),1-(row_len+row_space)*i},axisEnab($xaxstr)={col_len*j+col_space*j,(col_len)*(j+1)}
			// still need to work on x spacing
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i-1),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
			if(j==0)
				ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
			else
				ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j+col_space,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
			endif			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
		endfor
	endfor
end

//****************************************************************************************
function AppendWavesDifAxes(ywvstr,xwvstr,panel)
	string ywvstr // semi-colon delimited string of waves to graph
	string xwvstr //  semi-colon delimited string of waves to graph
	variable panel // which panel
	
	string yaxstr = "y" + num2istr(panel), xaxstr = "x" + num2istr(panel)
//	print ywvstr, xwvstr
	appendtograph/l=$yaxstr/b=$xaxstr $ywvstr vs $xwvstr
end

//****************************************************************************************
function AppendWavesDifAxesSubrange(ywvstr,xwvstr,panel,index)
	string ywvstr // semi-colon delimited string of waves to graph
	string xwvstr //  semi-colon delimited string of waves to graph
	variable panel // which panel
	variable index // the subrange index to graph
	
	string yaxstr = "y" + num2istr(panel), xaxstr = "x" + num2istr(panel)
	wave ywv = $ywvstr
	wave xwv = $xwvstr
	string tracename = ywvstr + "_" + num2istr(index)
	appendtograph/l=$yaxstr/b=$xaxstr ywv[index][]/TN=$tracename vs xwv
	ModifyGraph rgb($tracename)=(43520,43520,43520)
end