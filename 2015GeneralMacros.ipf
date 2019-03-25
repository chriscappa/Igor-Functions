#pragma rtGlobals=1		// Use modern global access method.
#pragma version=1.0	// to help with version control

#include <Axis Utilities>
#include<Remove Points>
#include<Concatenate Waves>

// *************************************************
//
//	2015GeneralMacros.ipf
//	Author:Donna Sueper
//	Modified by Ken Aikin
//    Revision Number: 1.00
//	Revision Date: Feb 19, 2015
//	Revision Notes: 


//	1.00	Started with the last version of 2013GeneralMacros. Added functions AddOrRemoveColon, Create2Dfrom1DFromTable and Create2Dfrom1DFromDataFolder from D.S.
//	20150219

// *************************************************

Menu "Macros"  // prevent the macros in this ipf from cluttering up the Macros menu.
End


// *************************************************

Menu  "2015General"

	"-"
	"Make Line On Graph"
	Submenu "Average"
		"AveragUsingStartStop"
		"AveragOnlyUsingStartStop"
		"AveragOnlyUzStartStopMinNonNans"
		"AveragXSecs"
		"AveragOnlyXSecs"
		"AveragOnlyXSecsMinNonNans"
		"Avg_WindDir_Xsec"
		"Calc_Avg_Alt_Profile"
	end
	Submenu "Create legends"
		"makeFzSizeLegend"
		"Color_legend_manual"
		"Color_legend_auto"
		"Size_legend_manual"
		"Size_legend_auto"
	end
	Submenu "Graph-Cursor"
		"Make_ones_between_cursors /1"
		"Make_nans_between_cursors /3"
		"Make_interp_between_cursors /5"
		"Make_avg_between_cursors /7"
		"Insert_new_avg /9"
		"Make_zeros_between_cursors /0"
	End
	Submenu "Kill All"
		"KillAllTables"
		"KillAllGraphs"
		"KillAllLayouts"
		"KillAllTablesGraphsLayoutsWaves"
		"KillAbsolutelyEverything"
	end
	Submenu "TimeExpansion"
		"Time_Expand"
		"Time Expand CreateWave"
		"Time Expand PrintMissing"
		"Time Expand PrintMissCreateWave"
		"Time Expand Dont Round Time"
		"Time Expand Dont Round CreateWave"
		"Time Expand Find Nearest"
		"Time Expand Fill"
	end
	"-"

End

// *************************************************
// Averages over segments of initial when key =1 and puts the segment averages in result.
// Numskip determines how many points after key=1 to ignore when averaging a segment.

Function Averag(result,initial,key,numskip)	
	Wave result,initial,key
	Variable numskip

	variable i=0, k=0, q=0, n=numpnts(initial)
	
	result = nan
	
	do
		if ((0.9<=key[i])*(key[i]<=1.1))
			k=i							//  i indicates the beginning of the mode, k  steps through the mode
			do							// skip first "numskip" seconds after mode begins (key goes to one)
				k+=1
			while ( (k< i + numskip)* (k<n) )
			q=k-1
			if ( (0.9<=key[k])*(key[k]<=1.1) )   	// make sure you are still in the mode
			//	q=k							// q marks the beginning of the averaging interval; k, the end
				do							// find k
					k+=1
				while ( (0.9<=key[k])*(key[k]<=1.1)*(k<n) ) 
				k-=1						// you will have gone one too many
				if (  (k-q)  >= 10)			//* make sure at least 10 values are in this interval
					WaveStats/q/r=[k, q] initial
					result[q+(k-q-1)/2]=v_avg
				endif
				i=k
			else					//* if not in measure mode after numskip-many points after the beginning of mode
				i=k
			endif	
		else
			i+=1
		endif
	while (i<n)

End


//*************************************************
//  This function averages values in TargetWave near transition points in the TransititionWave
// Input parameters:
// TargetWave is the wave to be averaged, typically a counts wave
// TransititionWave is typically a zero_key or cal-mode type of wave
// DestWave is the resulting wave - only this wave is changed by this function
// TransitionType is a string equal to one of the four values = "0to1", "1to0", "nantonum", or "numtonan" (all  capitalized/partially capilatized forms accepted)
// ForwardOrBackward  is a string equal to one of the two values = "forward" or "backward"  (all capitalized/partially capilatized forms accepted)
// Forward or backward refer to the direction one wishes to average in reference to the transition event.
// When "forward"  DeltaNum is the number of points to skip from the transition point in the forward direction.
// When "backward"  DeltaNum is the number of points to skip from the transition point in the backwards direction.
// NumAvgInterval in the number of points to average past the DeltaNum
// MinNumInAvg is the minimum number of non-nan values in the TargetWave  for which an average is deemed valid, and thus used.
//
// 	//To get a feel for what this function does, try executing the following lines of code: (decommentize first)
// make/o/n= 200 keywaveA, keywaveB, datawaveA, datawaveB, wave_avg1,  wave_avg2, wave_avg3,  wave_avg4, wave_avg5
// Edit datawaveA, keywaveA, wave_avg1,  wave_avg2, wave_avg3
// wave_avg1 = nan; wave_avg2 = nan; wave_avg3 = nan; wave_avg4 = nan; wave_avg5 = nan; 
// datawaveA = p
// keywaveA = SelectNumber( mod(p, 20) < 8, 0, 1)  // make a key wave with 8 1s and 12 zeros
//	 AveragAtTransition(datawaveA,keywaveA, wave_avg1, "0to1", "forward", 2, 2, 2)
//	 AveragAtTransition(datawaveA,keywaveA, wave_avg2, "0to1", "BAckward", 0, 2, 1)
//	 AveragAtTransition(datawaveA,keywaveA, wave_avg3, "1TO0", "forward", 3, 2, 1)
//
// 	//now with nans in the keywave, and what the heck, nans in the datawave
// keywaveB=keywaveA/keywaveA		// makes nans where keywaveA was zero
// datawaveB = selectnumber(mod(p, 3)!= 0, nan, datawaveA[p])		// makes nans at every 3rd point
// Edit datawaveB, keywaveB,  wave_avg4, wave_avg5
//
//  AveragAtTransition(datawaveB,keywaveB, wave_avg4, "NanToNum", "Forward", 3, 2, 1)
//  AveragAtTransition(datawaveB,keywaveB, wave_avg5, "numtonan", "backward", 2, 2, 2)
//
// 	// you can also decommentize the print statements in the DoAveragAtTransition for even more details
// // Users need to nan the DestWave before this function is called

Function AveragAtTransition(TargetWave,TransitionWave, DestWave, TransitionTypeStr, ForwardOrBackwardStr, DeltaNum, NumAvgInterval, MinNumInAvg)
	wave/d  TargetWave,TransitionWave, DestWave
	string  TransitionTypeStr, ForwardOrBackwardStr
	variable DeltaNum, NumAvgInterval, MinNumInAvg

	variable n							// a number of points in waves
	variable TransitionType    			// a numerical representation of the type
	variable ForwardOrBackward 		// a numerical representation of the direction
	variable transitionLeft, transitionRight		// will be passed to the DoAveragAtTransition function

	string TransitionTypeList = "0to1;1to0;nantonum;numtonan"		// make all lower case; user input get decapitalized
	string ForwardOrBackwardList = "forward;backward"  				// make all lower case; user input get decapitalized
	
	// some sanity checks
	
	n = numpnts(TargetWave)
	if (  (n != numpnts(TargetWave)) ||  n != numpnts(DestWave))
		abort "The number of points in all the waves must be the same - aborting"
	endif
	
	TransitionTypeStr = LowerStr(TransitionTypeStr)
	ForwardOrBackwardStr = LowerStr(ForwardOrBackwardStr)
	
	TransitionType = WhichListItem(TransitionTypeStr, TransitionTypeList, ";" )
	ForwardOrBackward = WhichListItem(ForwardOrBackwardStr, ForwardOrBackwardList, ";" )
	
	if (  NumAvgInterval < MinNumInAvg )
		abort "The value of NumAvgInterval needs to be >= MinNumInAvg  - aborting"
	endif

	if (  (TransitionType == -1) || (ForwardOrBackward== -1)  )
		abort "Please check the spelling of the transition type and the forward or backward parameter"
	endif

	
	// finished with sanity checks, get things set up
	
//	Users need to nan the wave themselves before this function is called
//	DestWave = nan  // don't do this because sometimes the same wave may need average results before and after a transition

	// we make a surrogate Transition wave for 2 reasons
	// 1.  The user may have generated the transition wave via some math, and some sneaky roundoff error could throw us off
	// 2.  The code is more compact if we use transition waves with 0s and 1s - it eliminates requent numtype() function calls.
	// If the user has a wave called MySurrogateTransitionWave already in the current directory,  too darn bad, it will be overwritten then killed.
	
	Duplicate/o TransitionWave MySurrogateTransitionWave
	MySurrogateTransitionWave  = round(TransitionWave)
	if (	TransitionType>=2) 	// we are dealing with nans - make a surrogate TransitionWave
		MySurrogateTransitionWave = selectnumber(numtype(MySurrogateTransitionWave[p])== 0, 0, 1)  // has 1 when a #, has a 0 when nan
	endif
	MySurrogateTransitionWave  = round(MySurrogateTransitionWave)

	switch(TransitionType)	// numeric switch
		case 0:			// 0 to 1
			transitionLeft=0
			transitionRight=1
			break
		case 1:			// 1 to 0
			transitionLeft=1
			transitionRight=0
			break
		case 2:			//nantonum
			transitionLeft=0 	// note that numtype(nan) == 2;  see MySurrogateTransitionWave above
			transitionRight=1 	// note that numtype(3) == 0;  see MySurrogateTransitionWave above
			break
		case 3:			// numtonan
			transitionLeft=1
			transitionRight=0
			break					
	endswitch
	
 	DoAveragAtTransition(TargetWave,MySurrogateTransitionWave, DestWave, ForwardOrBackward, DeltaNum, NumAvgInterval, MinNumInAvg, transitionLeft, transitionRight)
 	
	killwaves/Z MySurrogateTransitionWave

End

// *************************************************
// Same as AverageUsingStartStop except that no standard deviation or numpnts waves are created.
// For Really Really large data sets, is a bit faster than AveragUsingStartStop where we also create
// standard deviation and number of points waves.

Macro AveragOnlyUsingStartStop(timeWaveStr, wave2avgStr, waveAvgStr,startWaveStr, stopWaveStr)
	String timeWaveStr, wave2avgStr, waveAvgStr="my_avg",startWaveStr, stopWaveStr
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt waveAvgStr, "Name to give the averaged wave?"
	prompt startWaveStr, "Which wave has the start times?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt stopWaveStr, "Which wave has the stop times?",popup SortList(WaveList("*",";",""), ";", 16)

	if (numpnts($timeWaveStr)!=numpnts($wave2avgStr))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+ " - Aborting from AveragOnlyUsingStartStop"
	endif	
	
	if (numpnts($startWaveStr)!=numpnts($stopWaveStr))
		Abort " The start and stop waves must have the same number of points."+ startWaveStr+" "+stopWaveStr + " - Aborting from AveragOnlyUsingStartStop"
	endif	
	
	if (strlen(waveAvgStr)==0) 
		Abort "You need to insert values into the name prompts - Aborting from AveragOnlyUsingStartStop"
	endif
		
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=(Numpnts($startWaveStr)) $waveAvgStr
	
	DoAveragOnlyUsingStartStop($timeWaveStr, $wave2avgStr, $waveAvgStr, $startWaveStr, $stopWaveStr)
 
 End


// *************************************************
// Same as AveragOnlyUsingStartStop except that a minimum number
// of non-nan values must be present in the start-stop interval for an average to be calculated.

Macro AveragOnlyUzStartStopMinNonNans(timeWaveStr, wave2avgStr, waveAvgStr,startWaveStr, stopWaveStr, minNumNonNans)
	String timeWaveStr, wave2avgStr, waveAvgStr="my_avg",startWaveStr, stopWaveStr
	Variable minNumNonNans=1
	prompt timeWaveStr,"Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup WaveList("*",";","")
	prompt waveAvgStr, "Name to give the averaged wave?"
	prompt startWaveStr, "Which wave has the start times?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt stopWaveStr, "Which wave has the stop times?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt minNumNonNans, "Minimum number of non-nan values in each interval?"

	if (numpnts($timeWaveStr)!=numpnts($wave2avgStr))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragOnlyUzStartStopMinNonNans"
	endif	
	
	if (numpnts($startWaveStr)!=numpnts($stopWaveStr))
		Abort " The start and stop waves must have the same number of points."+ startWaveStr+" "+stopWaveStr+" - Aborting from AveragOnlyUzStartStopMinNonNans"
	endif	
	
	if (strlen(waveAvgStr)==0) 
		Abort "You need to insert values into the name prompts - Aborting from AveragOnlyUzStartStopMinNonNans"
	endif
		
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=(Numpnts($startWaveStr)) $waveAvgStr
	
	DoAveragOnlyStartStopMinNonNans($timeWaveStr, $wave2avgStr, $waveAvgStr, $startWaveStr, $stopWaveStr, minNumNonNans)
 
 End


// *************************************************
// Similar to AveragUsingStartStop only the start and stop times are defined by the number of seconds
// Waves start_xSecs and stop_xSecs are created and begin and end on "even" seconds
// For example if the timeWaveStr[0] = 12:14:08 and x = 60, 
// start_60sec[0] = 12:14:00, stop_60sec[0] = 12:14:59

Macro AveragOnlyXSecs(timeWaveStr, wave2avgStr, waveAvgStr, XSec)
	String timeWaveStr, wave2avgStr, waveAvgStr= "my_avg"
	Variable XSec=60
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt waveAvgStr,  "Name of the averaged wave?"
	prompt XSec, "Seconds to average over?"

	Silent 1

	if (numpnts($timeWaveStr)!=numpnts($wave2avgStr))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragOnlyXSecs"
	endif	
		
	if (strlen(waveAvgStr)==0)
		Abort "You need to insert values into the name prompts - Aborting from AveragOnlyXSecs"
	endif
		
	if (XSec<=0)
		Abort "XSec must be greater than zero."
	endif
	
	variable startTime =$timeWaveStr[0] - mod($timeWaveStr[0], XSec)
	variable stopTime = $timeWaveStr[numpnts($timeWaveStr)-1] -  mod($timeWaveStr[numpnts($timeWaveStr)-1], XSec ) + XSec
	
	Make/o/d/n=((stopTime-startTime)/XSec) $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	SetScale d 0,0,"dat", $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	$("start_"+num2str(XSec)+"sec") = startTime+p*XSec
	$("stop_"+num2str(XSec)+"sec") = startTime+(p+1)*XSec - 1
	
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=((stopTime-startTime)/XSec) $waveAvgStr
	
	DoAveragOnlyUsingStartStop($timeWaveStr, $wave2avgStr, $waveAvgStr, $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec"))

End


// *************************************************
// Creates averages based on number of seconds (default is 60 seconds) 

Macro AveragOnlyXSecsMinNonNans(timeWaveStr, wave2avgStr, waveAvgStr, XSec, minNumNonNans)
	String timeWaveStr, wave2avgStr, waveAvgStr= "my_avg"
	Variable XSec=60, minNumNonNans
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt waveAvgStr,  "Name of the averaged wave?"
	prompt XSec, "Seconds to average over?"
	prompt minNumNonNans, "Minimum number of non-nan values in each interval?"

	Silent 1

	if (minNumNonNans>XSec)
		Abort "Please choose values so that the number of seconds to average is greater than the minimum number of non-nan values  - Aborting from AveragOnlyXSecsMinNonNans"
	endif
	
	variable startTime = $timeWaveStr[0] - mod($timeWaveStr[0], XSec)
	variable stopTime = $timeWaveStr[numpnts($timeWaveStr)-1] -  mod($timeWaveStr[numpnts($timeWaveStr)-1], XSec ) + XSec
	
//	wave startwave =  $("start_"+num2str(XSec)+"sec")
	if(waveexists($("start_"+num2str(XSec)+"sec"))==0)
		Make/o/d/n=((stopTime-startTime)/XSec) $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
		SetScale d 0,0,"dat", $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	else
		redimension/n=((stoptime-starttime)/xsec)  $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	endif
	
	$("start_"+num2str(XSec)+"sec") = startTime+p*XSec
	$("stop_"+num2str(XSec)+"sec") = startTime+(p+1)*XSec - 1
	
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=((stopTime-startTime)/XSec) $waveAvgStr
	
	DoAveragOnlyStartStopMinNonNans($timeWaveStr, $wave2avgStr, $waveAvgStr, $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec"), minNumNonNans)

End


// *************************************************
// Averages data in wave2avgStr over the time interval indicated by 
// startWaveStr, stopWaveStr.
// It creates or overwrites waves with the names given in
// waveAvgStr,waveStdStr, and waveNpntsStr.
// It does not assume that data is in one second resolution. 
//  Assumes that all time waves are in increasing order.
// It does not round the time waves off to seconds.

Macro AveragUsingStartStop(timeWaveStr, wave2avgStr, waveAvgStr,waveStdStr, waveNpntsStr, startWaveStr, stopWaveStr)
	String timeWaveStr, wave2avgStr, waveAvgStr= "my_avg",waveStdStr= "my_std", waveNpntsStr= "my_num", startWaveStr, stopWaveStr
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt waveAvgStr,  "Name to give the averaged wave?"
	prompt waveStdStr,  "Name to give the std dev wave?"
	prompt waveNpntsStr,  "Name to give the number of points wave?"
	prompt startWaveStr, "Which wave has the start times?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt stopWaveStr, "Which wave has the stop times?",popup SortList(WaveList("*",";",""), ";", 16)
	
	if (numpnts($timeWaveStr)!=numpnts($wave2avgStr))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragUsingStartStop"
	endif	
	
	if (numpnts($startWaveStr)!=numpnts($stopWaveStr))
		Abort " The start and stop waves must have the same number of points."+ startWaveStr+" "+stopWaveStr+" - Aborting from AveragUsingStartStop"
	endif	
		
	if(    (strlen(waveAvgStr)==0) || (strlen(waveStdStr)==0) || (strlen(waveNpntsStr)==0)   )
		Abort "You need to insert values into the name prompts - Aborting from AveragUsingStartStop"
	endif	
	
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	waveStdStr = CleanUpName(waveStdStr, 0)
	waveNpntsStr = CleanUpName(waveNpntsStr, 0)
	Make/o/d/n=(numpnts($startWaveStr)) $waveAvgStr, $waveStdStr, $waveNpntsStr

	DoAveragUsingStartStop($timeWaveStr, $wave2avgStr, $waveAvgStr, $waveStdStr, $waveNpntsStr, $startWaveStr, $stopWaveStr)
 
 End


// *************************************************
// Similar to AveragUsingStartStop only the start and stop times are defined by the number of seconds
// Waves start_xSecs and stop_xSecs are created and begin and end on "even" seconds
// For example if the timeWaveStr[0] = 12:14:08 and x = 60, 
// start_60sec[0] = 12:14:00, stop_60sec[0] = 12:14:59

Macro AveragXSecs(timeWaveStr, wave2avgStr, waveAvgStr,waveStdStr, waveNpntsStr, XSec)
	String timeWaveStr, wave2avgStr, waveAvgStr= "my_avg",waveStdStr= "my_std", waveNpntsStr= "my_num"
	Variable XSec=60
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2avgStr, "Data to be averaged?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt waveAvgStr,  "Name of the averaged wave?"
	prompt waveStdStr,  "Name of the std dev wave?"
	prompt waveNpntsStr, "Name of the numpts wave?"
	prompt XSec, "Seconds to average over?"

	Silent 1

	if (numpnts($timeWaveStr)!=numpnts($wave2avgStr))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragXSecs"
	endif	
		
	if(    (strlen(waveAvgStr)==0) || (strlen(waveStdStr)==0) || (strlen(waveNpntsStr)==0)   )
		Abort "You need to insert values into the name prompts - Aborting from AveragXSecs"
	endif
	
	if (XSec<=0)
		Abort "XSec must be greater than zero."
	endif

	variable startTime =  $timeWaveStr[0] - mod( $timeWaveStr[0],XSec )
	variable stopTime =  $timeWaveStr[numpnts( $timeWaveStr)-1] -  mod( $timeWaveStr[numpnts( $timeWaveStr)-1],XSec ) + XSec
	
	Make/o/d/n=((stopTime-startTime)/XSec) $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	SetScale d 0,0,"dat",  $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	$("start_"+num2str(XSec)+"sec") = startTime+p*XSec
	$("stop_"+num2str(XSec)+"sec") = startTime+(p+1)*XSec - 1

	waveAvgStr = CleanUpName(waveAvgStr, 0)
	waveStdStr = CleanUpName(waveStdStr, 0)
	waveNpntsStr = CleanUpName(waveNpntsStr, 0)
	Make/o/d/n=((stopTime-startTime)/XSec) $waveAvgStr, $waveStdStr, $waveNpntsStr
	
	DoAveragUsingStartStop($timeWaveStr, $wave2avgStr, $waveAvgStr, $waveStdStr, $waveNpntsStr,  $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec"))

End


// *************************************************
//  Averages values in wave2avg at non-nan values in avg_key.
//  The width of one half the averaging interval is delta.
// If there is more than one non-nan value within delta of an earlier one,
// the function will skip over the second.

Function Averag_about_pt(waveAvg, wave2avg, avg_key, delta)
	Wave waveAvg, wave2avg, avg_key
	Variable delta

	if (numpnts(waveAvg)!= numpnts(wave2avg))
		Abort "The number of points in the two input waves of Averag_about_pt must be the same - Aborting from Averag_about_pt"
	endif
	
	variable n=numpnts(waveAvg), i=0
	
	i=FindFirst(avg_key, 0)
	
//	waveAvg=nan
	
	do
		if ( (i-delta >=0)&&(i+delta < n )  )
			WaveStats/q/r=[i-delta, i+delta] wave2avg
			waveAvg[i]=V_avg		
		endif
		i=FindFirst(avg_key, i+delta)
	while(i<n)

End


// *************************************************
//  Averages values in wave2avg according to where the state key  has successive values of 1.
// This functionality has been replaced in successive reduction schemes.  
// It is kept here for backwards compatibility.

Function Averag_about_state(waveAvg, wave2avg, state_key, meas_key, delta)
	Wave waveAvg, wave2avg, state_key,meas_key
	Variable delta

	variable n=numpnts(waveAvg), i=0
	
	//waveAvg=nan
	
	do
	
		do									// find the point where state_key = 1
			i+=1
		while ( (i<n)%&(state_key[i]==0) )
		
		do									// find the point where  meas_key = 1
			i-=1
		while ( (i>0)%&(meas_key[i]!=1) )
	
		if ( (i<(n-1)) %& (i-delta>=0) )
			WaveStats/q/r=[i-delta+1, i] wave2avg
			waveAvg[i-delta+1, i]=nan
			waveAvg[i+1-delta/2]=V_avg		
		endif
	
		do									// find the point where state_key = 0
			i+=1
		while ( (i<n)%&(state_key[i]==1) )
		
		do									// find the point where meas_key = 1
			i+=1
		while ( (i<n)%&(meas_key[i]!=1) )
		
		if ( (i+delta-1) < n )
			WaveStats/q/r=[i, i+delta-1] wave2avg
			waveAvg[i, i+delta-1]=nan
			waveAvg[i+delta/2]=V_avg		
		endif
		
	while(i<n)

End


//***************************************************
// Does the work of AveragAtTransition
Function DoAveragAtTransition(TargetWave,TransitionWave, DestWave, ForwardOrBackward, DeltaNum, NumAvgInterval, MinNumInAvg, transitionLeft, transitionRight)
	wave/d  TargetWave,TransitionWave, DestWave
	variable ForwardOrBackward, DeltaNum, NumAvgInterval, MinNumInAvg, transitionLeft, transitionRight

	variable x, n, startpt, stoppt, val, valnext
	variable transitionLength

	n = numpnts(TargetWave)

	for(x=0;  x<(n-1) ;  x+=1)	// assume a transition state at the end of a wave is invalid
		
		val = TransitionWave[x]
		valnext= TransitionWave[x+1]
			
		if ( (val == transitionLeft) && (valnext ==transitionRight)  )  	// this is a transition!
//			print x, x+1
			
			if (ForwardOrBackward == 0)	// FORWARD!
				startpt = x+1			// startpt is the first transitionRight that has a previous transitionLeft
				do
					x+=1
				while( (x< (n-1)) && (TransitionWave[x+1]==transitionRight) &&((x - (startpt+DeltaNum) + 1) < NumAvgInterval)  )
				
				stoppt = x  			//stoppt is the last 1 before a different value
				transitionLength = stoppt - (startpt+ DeltaNum) + 1
				
//				print startpt, stoppt, transitionLength		
				
				if ((transitionLength == MinNumInAvg) && (MinNumInAvg==1)) //Wavestats chokes when there is only one point in the data range
					DestWave[startpt]=TargetWave[startpt]
				else											 // we can average 

					if (transitionLength >=NumAvgInterval) // Note that transitionLength whould never be > NumAvgInterval, but I use >= instead of == just cuz
						Wavestats/q/r=[startpt+DeltaNum, stoppt] TargetWave
						if (V_npnts>=MinNumInAvg)
							DestWave[startpt+DeltaNum + transitionLength/2]=V_avg
						endif
					endif // transitionLength is long enough

				endif // transitionLength = 1 case
			
			else								// BACKWARD!	  beware, startpt > stoppt
				startpt = x  	// startpt is the first transitionLeft that has a following transitionRight
				do
					x-=1
				while( (x>= 0) && (TransitionWave[x]==transitionLeft) &&(startpt - (x+DeltaNum) < NumAvgInterval)  )
				
				stoppt = x +1 			//stoppt is the last transitionLeft before a different value
				transitionLength = (startpt- DeltaNum) - stoppt + 1

				//print startpt, stoppt, transitionLength		

				if ((transitionLength == MinNumInAvg) && (MinNumInAvg==1)) //Wavestats chokes when there is only one point in the data range
					DestWave[startpt]=TargetWave[startpt]
				else											 // we can average 

					if (transitionLength >=NumAvgInterval) // Note that transitionLength should never be > NumAvgInterval, but I use >= instead of == just cuz
						Wavestats/q/r=[startpt-DeltaNum, stoppt] TargetWave
						if (V_npnts>=MinNumInAvg)
							DestWave[stoppt + transitionLength/2]=V_avg		
//							print stoppt + transitionLength/2			
						endif					
					endif	// transitionLength is long enough
				
				endif  // transitionLength = 1 case
			
				x=startpt+1	//keep moving forward in the wave
			endif		// endif forward or backward
			
					
		else // we are not at a transition, move forward via the increment in the for loop

		endif
	
	endfor
End


// *************************************************
// Same as DoAveragOnlyUsingStartStop except that a minimum number of non-Nans in an interval must exist for a valid average

Function DoAveragOnlyStartStopMinNonNans(timeWave, wave2avg, waveAvg, startWave, stopWave, minNumNonNans)
	Wave timeWave, wave2avg, waveAvg,startWave, stopWave
	Variable minNumNonNans
	
	variable start=0, stop=0, i=0,numpts_start=numpnts(startWave)
	
	WaveStats/m=2/q timeWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (timeWave) + " has nans.  Cannot do average using start, stop - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
	WaveStats/m=2/q startWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (startWave) + " has nans.  Cannot do average using start, stop  - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
	WaveStats/m=2/q stopWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (stopWave) + " has nans.  Cannot do average using start, stop  - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
		
	waveAvg=nan
	
	start = 0
	stop = 0

	for(i=0;i<numpts_start;i+=1)	//bug. Was x<numpts_start. Changed to i<numpts_start	//ka 20070111
		start = FindFirstGE(timeWave, startWave[i], start)	
	  	stop = findRightMostLE(timeWave, stopWave[i], stop)	

		if ( (timeWave[start]>=startWave[i])&&(timeWave[stop]<=stopWave[i])  )			
			if ( (stop - start)>=1)
				WaveStats/m=2/q/r=[start,stop] wave2avg
				if (V_npnts>=minNumNonNans)				
					waveAvg[i] = V_avg
//					waveAvg[i]=calcAvg(start, stop, minNumNonNans, wave2avg)	//V_avg
				endif		
			else
				if (  ( (stop-start)==0)&&(minNumNonNans==1)	)		// There is only one point, start = stop.
					waveAvg[i]=wave2avg[start]
				endif
			endif		
		endif
		
	endfor

End

Function DoAveragStdStartStopMinNonNans(timeWave, wave2avg, waveAvg, waveStd, startWave, stopWave, minNumNonNans)
	Wave timeWave, wave2avg, waveAvg,waveStd, startWave, stopWave
	Variable minNumNonNans
	
	variable start=0, stop=0, i=0,numpts_start=numpnts(startWave)
	
	WaveStats/m=1/q timeWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (timeWave) + " has nans.  Cannot do average using start, stop - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
	WaveStats/m=1/q startWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (startWave) + " has nans.  Cannot do average using start, stop  - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
	WaveStats/m=1/q stopWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (stopWave) + " has nans.  Cannot do average using start, stop  - Aborting from DoAveragOnlyStartStopMinNonNans"
	endif
		
	waveAvg=nan
	waveStd=nan
	
	start = 0
	stop = 0

	for(i=0;i<numpts_start;i+=1)	//bug. Was x<numpts_start. Changed to i<numpts_start	//ka 20070111
		start = FindFirstGE(timeWave, startWave[i], start)	
	  	stop = findRightMostLE(timeWave, stopWave[i], stop)	

		if ( (timeWave[start]>=startWave[i])&&(timeWave[stop]<=stopWave[i])  )			
			if ( (stop - start)>=1)
				WaveStats/q/r=[start,stop] wave2avg
				if (V_npnts>=minNumNonNans)				
					waveAvg[i]=V_avg
					waveStd[i]=V_sdev
//					print waveAvg[i]
				endif		
			else
				if (  ( (stop-start)==0)&&(minNumNonNans==1)	)		// There is only one point, start = stop.
					waveAvg[i]=wave2avg[start]
				endif
			endif		
		endif
		
	endfor

End

function calcAvg(startRow, stopRow, minN, wname)	//pass in starting row, ending row, minimum N to calc a valid avg, and wave to average. Return avg.
	variable startRow, stopRow, minN	//if minN != 0, then consider its value
	wave wname
	
	variable i, n, theSum, avg
	
	n=0
	theSum=0
	avg=0
	
	for(i=startRow; i<=stopRow; i+=1)
		if (numtype(wname[i])==0)	//if its a valid number
			theSum+= wname[i]		//sum numbers
			n+=1						//# of values going into avg
		endif
	endfor
	
	if (minN > 0)		//if minN is specified
		if (n >= minN)	//then do the avg
			avg = theSum/n
		else				//not enough valid points. No avg
			avg = NaN
		endif
	else					//don't care about minimum number of points to calc avg
		if (n>0)		//if there was at least one point, calc avg
			avg = theSum/n
		else				//no data to avg
			avg = NaN
		endif
	endif
	
	return avg
end function

// *************************************************
// Same as DoAveragUsingStartStop except the std dev and num pnts waves are NOT created.

Function DoAveragOnlyUsingStartStop(timeWave, wave2avg, waveAvg, startWave, stopWave)
	Wave timeWave, wave2avg, waveAvg,startWave, stopWave

	variable start=0, stop=0, i=0,numpts_start=numpnts(startWave)
	
	WaveStats/m=2/q timeWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (timeWave) + " has nans.  Cannot do average using start, stop - Aborting from DoAveragOnlyUsingStartStop"
	endif
	WaveStats/m=2/q startWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (startWave) + " has nans.  Cannot do average using start, stop - Aborting from DoAveragOnlyUsingStartStop"
	endif
	WaveStats/m=2/q stopWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (stopWave) + " has nans.  Cannot do average using start, stop - Aborting from DoAveragOnlyUsingStartStop"
	endif
		
	waveAvg=nan
		
	start = 0
	stop = 0

	for(i=0;i<numpts_start;i+=1)
		start = FindFirstGE(timeWave, startWave[i], start)	
	  	stop = findRightMostLE(timeWave, stopWave[i], stop)		
		
		if ( (timeWave[start]>=startWave[i])&&(timeWave[stop]<=stopWave[i])  )			
			if ( (stop - start)>=1)
				WaveStats/m=2/q/r=[start,stop] wave2avg
				waveAvg[i]=V_avg		
			else
				if ( (stop - start)==0)		// There is only one point.
					waveAvg[i]=wave2avg[start]
				endif
			endif		
		endif
		
	endfor

End


// *************************************************
// Calculates the average, standard deviation and number of non-nan points between time intervals start and stop times waves

Function DoAveragUsingStartStop(timeWave, wave2avg, waveAvg, wave_std, wave_numpnts, startWave, stopWave)
	Wave timeWave, wave2avg, waveAvg, wave_std, wave_numpnts,startWave, stopWave

	variable start, stop, n=numpnts(timeWave)	// start and stop range over the values 0 through numpnts(timeWave)
	variable i=0, numpts_start=numpnts(startWave) // i ranges over the values 0 through numpnts(startWave)
	
	WaveStats/q timeWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (timeWave) + " has nans - Aborting from DoAveragUsingStartStop"
	endif
	WaveStats/q startWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (startWave) + " has nans - Aborting from DoAveragUsingStartStop"
	endif
	WaveStats/q stopWave
	if (V_numNans !=0)
		Abort "The timewave "+NameofWave (stopWave) + " has nans - Aborting from DoAveragUsingStartStop"
	endif
		
	waveAvg=nan
	wave_std=nan
	wave_numpnts=nan

	start = 0
	stop = 0

	for(i=0;i<numpts_start;i+=1)
		start = FindFirstGE(timeWave, startWave[i], start)	
	  	stop = findRightMostLE(timeWave, stopWave[i], stop)		

		if ( (timeWave[start]>=startWave[i])&&(timeWave[stop]<=stopWave[i])  )			
			if ( (stop - start)>=1)
				WaveStats/q/r=[start,stop] wave2avg
				waveAvg[i]=V_avg		
				wave_std[i]=V_sdev
				wave_numpnts[i]=V_npnts
			else
				if ( (stop - start)==0)
					waveAvg[i]=wave2avg[start]
					wave_std[i]=0
					wave_numpnts[i]=1
				endif
			endif		
		endif
	endfor
			
End


// *************************************************
// Replaces the replaceNum values in a wave with Nan values.

Function Change2Nans(w,replaceNum)
	Wave w
	Variable replaceNum

// be aware of roundoff errors!!
// w = selectnumber(w[p]==replaceNum, w[p], nan) // this is not matrix aware

	variable i=0, n=numpnts(w)
	
	for(i=0; i<n;i+=1)
		if (w[i] ==replaceNum)
			w[i] =nan
		endif	
	endfor
	
End


// *************************************************
// Replaces the replaceNum values in a wave with Nan values.

Function Change2NansRound(w,replaceNum)
	Wave w
	Variable replaceNum

// w = selectnumber(w[p]==replaceNum, w[p], nan) // this is not matrix aware

	variable i=0, n=numpnts(w)
	
	for(i=0; i<n;i+=1)
		if (round(w[i]) == replaceNum)
			w[i] =nan
		endif	
	endfor
	
End


// *************************************************
// Replaces the replaceNum values for all the waves in the wavelist with Nan values.

Function Change2NansWList(wList,replaceNum)
	String wList
	Variable replaceNum

// be aware of roundoff errors!!

	variable i=0, n=ItemsInList(wList)
	string oneWave
	
	for(i=0;i<n;i+=1)
		oneWave = StringFromList(i, wList,";")
		if (exists(oneWave)==1)
			WAVE w = $oneWave
			Change2Nans(w,replaceNum)
		else
			Print "// The following is not a wave: " + oneWave+ "; printing from Change2NansWList"
		endif
	endfor

End


//*************************************
//  Replaces everything in the wave with a Nan or a number; changes only +-infs.

Function ChangeInfs2Nans(w)
	Wave w

	variable i, y,n=numpnts(w)
	
	for (i=0;i<n;i+=1)
		y = w[i]
		w[i] = selectnumber(numtype(y)==0, nan, y)		// don't use p for the index - make it general for 2 dimensional waves
	endfor
		
End


// *************************************************
// Replaces the NaN values in a wave with replaceNum.

Function ChangeNans(w,replaceNum)
	Wave w
	Variable replaceNum

// be aware of roundoff errors!!
//	w = selectnumber(numtype(w[p])==2, w[p], replaceNum)  // this is not matrii aware

	variable i=0, n=numpnts(w)
	
	for(i=0; i<n;i+=1)
		if (numtype(w[i]) ==2)
			w[i] = replaceNum
		endif	
	endfor
	
End


// *************************************************
// Replaces the NaN values for all the waves in the wavelist with the replaceNum.

Function ChangeNansWList(wList,replaceNum)
	String wList
	Variable replaceNum

// be aware of roundoff errors!!

	variable i=0, n=ItemsInList(wList)
	string oneWave
	
	for(i=0;i<n;i+=1)
		oneWave = StringFromList(i, wList,";")
		if (exists(oneWave)==1)
			wave w = $oneWave
			ChangeNans(w,replaceNum)
		else
			Print "// The following is not a wave: " + oneWave+ "; printing from Change2NansWList"
		endif
	endfor

End

// *************************************************
//  Duplicates data wave, only at points where state_key = 1 and 
// after stateBeginDelta number of points after the state has begun and 
// stateEndDelta number of points before the state ends.
// Can be used for example, to extract NO_counts values from the NO_zero_k state.
// Key state can have 0s or nans elsewhere, it just looks for 1s.
// Call this using something like ExtractDataFromState(NO_counts, NO_zero_k, "NO_zero_counts", 5, 2)
// Newly created wave is named NO_zero_counts.

Function ExtractDataFromState(dataWave, keyWave, newWave, stateBeginDelta, stateEndDelta)
	Wave dataWave, keyWave
	String newWave
	Variable stateBeginDelta, stateEndDelta

	if (numpnts(dataWave)!=numpnts(keyWave) )
		Abort "The number of points for the two input waves in ExtractDataFromState must be the same - Aborting from ExtractDataFromState"
	endif
	
	variable n=numpnts(dataWave), x=0, StateBeginPt=0, StateEndPt=0

	Make/o/d/n=(n) $(newWave)
	WAVE stateWave =$(newWave)
	stateWave = nan

	do

		do									// find the point where keyWave = 1
			StateBeginPt+=1
		while ( (StateBeginPt<n)&&(keyWave[StateBeginPt]!=1) )
			
		StateEndPt = StateBeginPt
		do									// find the point where  keyWave != 1
			StateEndPt+=1
		while ( (StateEndPt<n)&&(keyWave[StateEndPt]==1) )
		
		if (  (StateEndPt==n)&&(keyWave[n-1]==1)  )		// check to see if we ended in a 1 state
			StateEndPt = n-1
		endif

		if (   ( (StateBeginPt+stateBeginDelta)<StateEndPt) && (StateEndPt<n) && ( (StateBeginPt+stateBeginDelta) < (StateEndPt-stateEndDelta) )  )
			stateWave[StateBeginPt+stateBeginDelta, StateEndPt-stateEndDelta] = dataWave[p]	
		endif

		StateBeginPt = StateEndPt +1
	
	while(StateBeginPt<n)

End


// *************************************************
//  Puts numOfNans-many nans in wave2BeNaned whenever keyWave is 1.
//  This one is old legacy code for ignoring values when a state changed

Function NAN_according_to_key(wave2BeNaned, keyWave, numOfNans)	
	Wave wave2BeNaned, keyWave
	Variable numOfNans

	NanAccordingToVal(wave2BeNaned, keyWave, numOfNans, 1)

End


// *************************************************
//  Puts numOfNans-many nans in wave2BeNaned whenever keyWave is val.

Function NanAccordingToVal(wave2BeNaned, keyWave, numOfNans, val)	
	Wave wave2BeNaned, keyWave
	Variable numOfNans, val

	variable i=0, n=numpnts(wave2BeNaned)
	
	do			// if keyWave ~= val, assign nan to subsequent numOfNans seconds
		if(val>0)
			if(  (keyWave[i] >= val*0.99)&&(keyWave[i] <= val*1.01)  )
				wave2BeNaned[i, i+numOfNans]=nan
			endif
		else		// val is negative
			if(  (keyWave[i] <= val*0.99)&&(keyWave[i] >= val*1.01)  )
				wave2BeNaned[i, i+numOfNans]=nan
			endif
		endif	// note that if val is 0, above checks for equality, not within 1%
		
		i+=1
	while (i<n)

End


// *************************************************
// Returns the first non-nan value after and including beginPt in a wave.
// If all values after beginPt are nans, Findfirst will return number of points.

Function FindFirst(w, beginPt)
	Wave w
	Variable beginPt

	variable i=beginPt,  n=numpnts(w)

	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif
		
	if (numtype(w[i])== 0 )	// is w[beginPt] a number?
		return i	
	endif

	for (i=beginPt; i<n; i+=1)
		if (numtype(w[i])==0)
			break
		endif
	endfor
	
	return i

End
  

// *************************************************
// Returns the first point greater than or equal to value in the wave w.  
// It begins searching at beginPt and steps through the wave in INCREASING point number.
// It returns beginPt if w[beginPt] >=value.
//  It returns n if  all the points in w are less than value
//  It returns -1 if  beginPt < 0
// Typically used with waves that increase with increasing point number, like timewaves.

Function FindFirstGE(w, value, beginPt)
	Wave w
	Variable value, beginPt

	variable i=beginPt,  n = numpnts(w)
	
	if (i >= n)
		return n
	endif
	
	if (i < 0)
		return -1
	endif

	if (w[i] >= value)
		return i
	endif

	for(i=beginPt; i<n; i+=1)
		if (w[i]>=value)
			break
		endif
	endfor

	return i

End


// *************************************************
// Returns the first point less than or equal to value in the wave w. 
//  It begin searching at beginPt and steps through the wave in INCREASING point number.
// It returns beginPt if w[beginPt] <=value.
//  It returns n if  all the points are greater than value
//  It returns -1 if  beginPt < 0
// It returns i-1 if w[i] < value but  if w[i-1] > value.
// Typically used with waves that decrease with increasing point number.
// Rarely used when w is a timewave, because time usually runs forward.

Function FindFirstLE(w, value, beginPt)
	Wave w
	Variable value, beginPt

	variable i=beginPt,  n = numpnts(w)
	
	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (w[i] <= value)
		return i
	endif

	for(i=beginPt; i<n; i+=1)
		if (w[i]<=value)
			break
		endif
	endfor

	return i
	
End


// *************************************************
// Returns the first nan value after and including beginPt in a wave.
// If all values after beginPt are nans, Findfirst will return number of points.

Function FindFirstNan(w, beginPt)
	Wave w
	Variable beginPt

	variable i=beginPt, n=numpnts(w)

	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (numtype(w[i])==2)
		return i
	endif
	
	for(i=beginPt; i<n; i+=1)
		if (numtype(w[i])==2)
			break
		endif
	endfor
	
	return i

End


// *************************************************
// Returns the last non-nan value before (or including) endPt in a wave.
// If the entire wave is nans, FindLast will return -1.
// It steps through the wave backwards in decreasing point number.

Function FindLast(w, endPt)
	Wave w
	Variable endPt

	variable i=endPt, n=numpnts(w)

	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (numtype(w[i])==0)// is the w[endPt] a number?
		return i		
	endif

	for(i=endPt; i>=0; i-=1)
		if (numtype(w[i])==0)
			break
		endif
	endfor
	
	return i

End


// *************************************************
// Returns the last point greater than or equal to value in  wave w.  
// It begin searching at beginPt and steps through the wave in DECREASING point number
// It returns beginPt if w[beginPt] >= value but  w[beginPt-1] < value.
// Typically used with waves that increase with increasing point number, like timewaves.
// It steps through the wave backwards, in decreasing point number.

Function FindLastGE(w, value, beginPt)
	Wave w
	Variable value, beginPt

	variable i=beginPt, n= numpnts(w)
	
	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (w[i] >= value)// w[beginPt] is less than or equal to value
		return i
	endif

	 // w[i] <value
	for(i=beginPt; i>=0; i-=1)
		if (w[i]>=value)
			break
		endif
	endfor

	return i

End


// *************************************************
// Returns the last point less than or equal to value in wave w. 
// It returns beginPt if w[beginPt] <= value but  if w[beginPt-1] > value.
// Rarely used when w is a timewave, because time usually runs forward.
// It steps through the wave backwards; in decreasing point number.

Function FindLastLE(w, value, beginPt)
	Wave w
	Variable value, beginPt

	variable i=beginPt, n= numpnts(w)
	
	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (w[i] <= value) // w[beginPt] is less than or equal to value
		return i
	endif
	
	for(i=beginPt; i>=0; i-=1)
		if (w[i]<=value)
			break
		endif
	endfor

	return i

End


// *************************************************
// Returns the last nan value before and including beginPt in a wave.
// If all values after beginPt are non-nans, FindLastNan will return -1.
// It steps through the wave backwards; in decreasing point number.

Function FindLastNan(w, endPt)
	Wave w
	Variable endPt

	variable i=endPt, n=numpnts(w)

	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (numtype(w[i])==2)
		return i
	endif

	for(i=endPt; i>=0; i-=1)
		if (numtype(w[i])==2)
			break
		endif
	endfor

	return i

End


// *************************************************
// Returns the right-most point less than or equal to value 
// in the wave w.  It begin searching at beginPt and continues in increasing point number.
// It returns beginPt if w[beginPt] > value.
// This function assumes that w is sorted in increasing order.

Function FindRightMostLE(w, value, beginPt)
	Wave w
	Variable value, beginPt

	variable i=beginPt, n= numpnts(w)
	
	if (i >= n)
		return n
	endif

	if (i < 0)
		return -1
	endif

	if (w[i] > value)
		return beginPt	//=i
	endif
	
	if (w[n-1] <= value)	// the last point in  w is less than or equal to value
		return (n-1)
	endif
	
	do
		if (w[i+1] <=value)
			i+=1
		endif
	while(   (i<(n-1)) && (w[i] <=value) && ( w[i+1] <= value )  )
		
	return i
	
End


// *************************************************
// Duplicates first non-nan ( after startpoint ) to the first point in the wave 
// and the last non-nan (before endPt) to the end of the wave
// Does no error checking

Function InsertBeginEndPoints(w, startpoint, endPt)
	Wave w
	Variable startpoint, endPt

	w[startpoint]=w[FindFirst(w,(startpoint+1))]
	w[endPt]=w[FindLast(w,(endPt-1))]

End


// *************************************************
// Performs a linear interpolation on wave originalWave.
// The interpolation built into igor can change wave scaling, so be careful if you use that.

Function InterpolateWave(resultWave, originalWave)		
	Wave resultWave, originalWave

	variable i, y, nOrig, nResult
	
	string cmd, nOrigStr
	string waveInterpolatedStr, wave2interpStr, iStr
	
	nOrig=numpnts(originalWave)
	nResult = numpnts(resultWave)
	
	if  (nOrig != nResult)
		Abort "The number of points in both input waves should be the same - Aborting from InterpolateWave"
	endif
	
	
	iStr= UniqueName("i", 1, 0)				// get wave names that aren't being used
	waveInterpolatedStr = UniqueName("w_i", 1, 0)		// get wave names that aren't being used
	wave2interpStr = UniqueName("w_2i", 1, 0)

	Make/d/o/n=(nOrig) $iStr, $waveInterpolatedStr,   $wave2interpStr
	WAVE iwave = $iStr
	WAVE w_i = $waveInterpolatedStr
	WAVE w_2i = $wave2interpStr
	iwave = p
	w_2i = originalWave
	WaveStats /Q w_2i
	if (V_numNans!=nOrig)
		i=FindFirst(w_2i,0)
		y=FindLast(w_2i,nOrig-1)
		w_2i[0]=w_2i[i]
		w_2i[nOrig-1]=w_2i[y]

		sprintf nOrigStr, "%d", nOrig	//do this instead of num2str. Gets all the sig figs,so that cmd operates correctly. Fixed bug 3/9/07.
		cmd="interpolate /T=1/N="+nOrigStr+"/y="+waveInterpolatedStr+", "+wave2interpStr + " /X="+iStr
//		print cmd
		Execute cmd

		resultWave=w_i
	else
		Print "// Interpolation is a problem for wave "+wave2interpStr
		resultWave = nan
	endif
	
	Killwaves iwave, w_i,w_2i

End


// *************************************************
// Duplicates the last non-nan value in a wave  to the last point in a wave
// if there is no last non nan value the function returns 0
// Useful for waves which indicatee state (zero, measure) but may not have a reading, non-nan value every second
// All nan values preceding the first nan value remain nan.

Function RepeatLastNonNan(w)
	Wave w

	variable i=FindFirst(w, 0), n=numpnts(w)
	variable temp = w[i]
	
	if (i == n)
		return 0		// the entire wave is nans
	endif
	
	do
		if(numtype(w[i])==2)
			w[i]=temp
		else
			temp = w[i]
		endif 	
		i+=1
	while(i<n)

	return 1
	
End


// *************************************************
//  Duplicate the columns of matrix w1 to matrix w2.  w2 will have the 
// number of columns equal to the sum of columns w1 and w2.

Function ConcatenateMatrixCols(w1, w2)
	String w1, w2
	
	variable numPoints1, numPoints2
	variable  w1numRows, w1numCols,  w2numRows, w2numCols

	string wInfo
	variable WType
	
	if (Exists(w1) == 0)
		Duplicate/R=[] $w2, $w1
	else
		wInfo=WaveInfo($w2, 0)
		wType=NumberByKey("NUMTYPE", wInfo)
		if (wType)				// Numeric wave
			WAVE/C/D wave1 = $(w1)
			WAVE/C/D wave2 = $(w2)
			w1numRows = DimSize(wave1, 0)
			w1numCols = DimSize(wave1, 1)
			w2numRows = DimSize(wave2, 0)
			w2numCols = DimSize(wave2, 1)	
			Redimension/N=(w1numRows, w1numCols + w2numCols) $w1
	
//			print w1numRows, w1numCols,  w2numRows, w2numCols
			if (w1numRows != w2numRows)
				Abort "the number of rows must be the same to concatenate column-wise -  Aborting from ConcatenateMatrixCols"
			endif
		
			wave1[, ][w1numCols,] = wave2[p][q-w1numCols]
		else						// Text wave
			WAVE/T twave1 = $(w1)
			WAVE/T twave2 = $(w2)
			w1numRows = DimSize(twave1, 0)
			w1numCols = DimSize(twave1, 1)
			w2numRows = DimSize(twave2, 0)
			w2numCols = DimSize(twave2, 1)		
			if (w1numRows != w2numRows)
				Abort "the number of rows must be the same to concatenate column-wise. -  Aborting from ConcatenateMatrixCols"
			endif
			twave1[w1numRows, ][w1numCols,] = twave2[p][q-w1numCols]
		endif
	endif
End


// *************************************************
// Transposes a table with 1 dimesional waves; creates new waves with the transposed values.
// Will not work for waves in different data folders.

Function TransposeTable()

	if (WinType("")!=2)
		Abort "The target window must be a table to use the transposeTable function. - Aborting from TransposeTable"
	endif
		
	if (strlen(StringByKey("COLUMNNAME",TableInfo("", 0) ))==0)
		Abort "Please call this function on a non-empty window. - Aborting from TransposeTable"
	endif
	
	variable colIndex, rowIndex, colInTable
	string tableWaveStr,tableWaveList
	variable WaveNumType, WaveNumPts, thisWaveNumType, thisWaveNumPts
	
	colInTable = NumberByKey("COLUMNS",TableInfo("", -2) ) - 1 // subtract point column
	WaveNumPts = NumberByKey("ROWS",TableInfo("", -2) )
	
	tableWaveStr = StringByKey("COLUMNNAME",TableInfo("", 0) )
	tableWaveStr =tableWaveStr[0, strlen(tableWaveStr)-3]
	tableWaveList= tableWaveStr+";"
	
	WaveNumType = NumberByKey( "NUMTYPE", waveinfo($tableWaveStr,0) )
	
	for (colIndex=1;colIndex<colInTable;colIndex+=1 )
		tableWaveStr = StringByKey("COLUMNNAME",TableInfo("", colIndex) )
		tableWaveStr =tableWaveStr[0, strlen(tableWaveStr)-3]
		tableWaveList+= tableWaveStr+";"
		
		thisWaveNumType = NumberByKey( "NUMTYPE", waveinfo($tableWaveStr,0) )
		thisWaveNumPts = numpnts($tableWaveStr)
		
		if (WaveNumPts != thisWaveNumPts)
			Abort "The number of points in all the waves must be the same. "+tableWaveStr+" Aborting from TransposeTable"
		endif
		if (   ( (WaveNumType==0)%^(thisWaveNumType==0) )    )
			Abort "You cannot use both text and numberic waves for this operation. - Aborting from TransposeTable"
		endif
	
	endfor
	
	string NewWaveNameStr, textWaveNameStr
	prompt NewWaveNameStr, "Enter the base name of the new waves you would like to create."
	prompt textWaveNameStr, "Alternatively, select a text wave with the wave names.", popup,  " ;"+TextWaveList("*",";","")
	doPrompt "Enter base name or select a text wave with the wave names" NewWaveNameStr, textWaveNameStr
	
	if (V_flag == -1)
		Abort "The user has canceled the operation. - Aborting from TransposeTable"
	endif
	
	if (   (strlen(NewWaveNameStr)==0)&&(cmpstr(textWaveNameStr, " ")==0)  )
		Abort "Please try again with a new wave name or a selected text wave. - Aborting from TransposeTable"
	endif
	
	if(strlen(NewWaveNameStr)>0)
		Make/o/t/n=(WaveNumPts) TableWaveNames
		TableWaveNames = NewWaveNameStr+num2str(p)
	else
		if (numpnts($textWaveNameStr)!=WaveNumPts)
			Abort "Please select a text wave with the same number of points as the waves in the table. - Aborting from TransposeTable"
		endif
		Duplicate/t/o $textWaveNameStr TableWaveNames
	endif
	
	if (WaveNumType==0)
		make/t/o/n=(WaveNumPts, colInTable) NewTextMatrixToTranspose
		for(colIndex=0;colIndex<colInTable ;colIndex+=1)	
			tableWaveStr = StringFromList( colIndex,tableWaveList, ";")
			Wave/t tw = $tableWaveStr
			NewTextMatrixToTranspose[][colIndex] = tw[p]
		endfor		
		MatrixTranspose NewTextMatrixToTranspose
		for(rowIndex=0;rowIndex<WaveNumPts ;rowIndex+=1)	
			make/t/o/n=(colInTable) $(TableWaveNames[rowIndex])
			wave/t tw = $(TableWaveNames[rowIndex])
			tw = NewTextMatrixToTranspose[p][rowIndex]
		endfor
		Killwaves NewTextMatrixToTranspose
	else
		Make/o/d/n=(WaveNumPts, colInTable) NewNumMatrixToTranspose
		for(colIndex=0;colIndex<colInTable ;colIndex+=1)	
			tableWaveStr = StringFromList(colIndex, tableWaveList, ";")
			Wave nw = $(tableWaveStr)
			NewNumMatrixToTranspose[][colIndex] = nw[p]
		endfor		
		MatrixTranspose NewNumMatrixToTranspose
		for(rowIndex=0;rowIndex<WaveNumPts ;rowIndex+=1)	
			Make/o/d/n=(colInTable) $(TableWaveNames[rowIndex])
			wave nw = $(TableWaveNames[rowIndex])
			nw = NewNumMatrixToTranspose[p][rowIndex]
		endfor
		Killwaves NewNumMatrixToTranspose
	endif

	Edit 
	for(rowIndex=0;rowIndex<WaveNumPts ;rowIndex+=1)	
		AppendToTable $TableWaveNames[rowindex]; DelayUpdate
	endfor
	
	Killwaves TableWaveNames

End


// *************************************************
// Returns the name of the text wave created by converting the timewave named  in timeWavestr into a string containing yymmdd....

Function/S MakeYMDHMS(timeWaveStr,textTimeWaveStr)
	String timeWaveStr, textTimeWaveStr

	WAVE timeWave = $timeWaveStr
	WAVE/T textTimeWave = $textTimeWaveStr
	
	if(WaveExists(timewave)==0)
		Abort "There is no wave named "+timeWaveStr+" - Aborting from MakeYMDHMS"
	endif
	
	if(WaveExists(textTimeWave)==0 )
		textTimeWaveStr = CleanupName(textTimeWaveStr, 0)
		Make/t/o/n=(numpnts(timewave))  $textTimeWaveStr
		WAVE/t textTimeWave = $textTimeWaveStr
	else 
		if (numpnts(timewave) != numpnts(textTimeWave) )
			Abort "The number of points must be the same for MakeYYMDHMS to work - Aborting from MakeYMDHMS"
		endif
	endif
	
	textTimeWave = Time2yymmddhhmmss(timeWave[p])
	
	return textTimeWaveStr
	
End


// *************************************************
// Returns the name of the text wave created by converting the timewave named  
// in timeWavestr into a string containing yyyymmddhhmmss.

Function/S MakeYYMDHMS(timeWaveStr,textTimeWaveStr)
	String timeWaveStr, textTimeWaveStr

	WAVE timeWave = $timeWaveStr
	WAVE/t textTimeWave = $textTimeWaveStr
	
	if(WaveExists(timewave)==0)
		Abort "There is no wave named "+timeWaveStr+"- Aborting from MakeYYMDHMS"
	endif
	
	if(WaveExists(textTimeWave)==0 )
		textTimeWaveStr = CleanupName(textTimeWaveStr, 0)
		Make/t/o/n=(numpnts(timewave))  $textTimeWaveStr
		WAVE/t textTimeWave = $textTimeWaveStr
	else 
		if (numpnts(timewave) != numpnts(textTimeWave) )
			Abort "The number of points must be the same for MakeYYMDHMS to work - Aborting from MakeYYMDHMS"
		endif
	endif
	
	textTimeWave = Time2yyyymmddhhmmss(timeWave[p])

	return textTimeWaveStr

End


// *************************************************
// Returns a string with a 2 char year, 2 char month, etc. from timeVar

Function/s Time2yymmddhhmmss(timeVar)
	Variable TimeVar

	string returnStr
	string dateStr, timeStr, yy,mm,dd,hh,mimi,ss
	variable slashPos1, slashPos2
	
	dateStr=Secs2Date(timeVar,-1)  // example: 18/04/2004"  dd/mm/yyyy

	if (strsearch(dateStr, "(", 0)>0 )		// normal
		slashPos1=strsearch(dateStr, "/", 0)
		if(slashPos1==1)
			dd = "0"+dateStr[0]
		else
			dd =dateStr[0,1]
		endif
		
		slashPos2=strsearch(dateStr, "/", slashPos1+1)
		
		if(  (slashPos2 - slashPos1)==2)
			mm =  "0"+dateStr[slashPos2-1]
		else
			mm =dateStr[slashPos1+1,slashPos2-1]
		endif
		yy=dateStr[strlen(dateStr)-6, strlen(dateStr)-5]	//  count from the end and ignore the last 4 characeters, space, (, x, )
	
	else			// some windows machines aren't cooperative
	
		slashPos1=strsearch(dateStr, "/", 0)
		if(slashPos1==1)
			mm = "0"+dateStr[0]
		else
			mm =dateStr[0,1]
		endif
		
		slashPos2=strsearch(dateStr, "/", slashPos1+1)
		
		if(  (slashPos2 - slashPos1)==2)
			dd =  "0"+dateStr[slashPos2-1]
		else
			dd =dateStr[slashPos1+1,slashPos2-1]
		endif
			
		yy=dateStr[slashPos2+1, strlen(dateStr) -1]

	endif
	
	timeStr=Secs2Time(timeVar,3)
	hh=timeStr[0, 1]
	mimi=timeStr[3,4]
	ss=timeStr[6,7]
		
	returnStr = yy+mm+dd+hh+mimi+ss
		
	return returnStr

End


// *************************************************
// Returns a string with a 4 char year, 2 char month, etc. from timeVar

Function/S Time2yyyymmddhhmmss(timeVar)
	Variable timeVar

	string returnStr
	string dateStr, timeStr, yyyy,mm,dd,hh,mimi,ss
	variable slashPos1, slashPos2
	
	dateStr=Secs2Date(timeVar,-1) 		//	dateStr =  "18/04/2004" 

	if (strsearch(dateStr, "(", 0)>0 )		// normal
		slashPos1=strsearch(dateStr, "/", 0)
		if(slashPos1==1)
			dd = "0"+dateStr[0]
		else
			dd =dateStr[0,1]
		endif
		
		slashPos2=strsearch(dateStr, "/", slashPos1+1)
		
		if(  (slashPos2 - slashPos1)==2)
			mm =  "0"+dateStr[slashPos2-1]
		else
			mm =dateStr[slashPos1+1,slashPos2-1]
		endif
		yyyy=dateStr[strlen(dateStr)-8, strlen(dateStr)-5]	//  count from the end and ignore the last 4 characeters, space, (, x, )
	
	else	 	// some windows machines aren't cooperative
	
		slashPos1=strsearch(dateStr, "/", 0)
		if(slashPos1==1)
			mm = "0"+dateStr[0]
		else
			mm =dateStr[0,1]
		endif
		
		slashPos2=strsearch(dateStr, "/", slashPos1+1)
		
		if(  (slashPos2 - slashPos1)==2)
			dd =  "0"+dateStr[slashPos2-1]
		else
			dd =dateStr[slashPos1+1,slashPos2-1]
		endif
			
		yyyy=dateStr[slashPos2+1, strlen(dateStr) -1]
		 if (cmpstr(yyyy[0], "0")==0)		// won't work in 2010+
		 	yyyy = "20"+yyyy
		 else
		 	yyyy = "19"+yyyy
		 endif

	endif
		
	timeStr=Secs2Time(timeVar,3)
	hh=timeStr[0, 1]
	mimi=timeStr[3,4]
	ss=timeStr[6,7]
		
	returnStr = yyyy+mm+dd+hh+mimi+ss
		
	return returnStr

End


// *************************************************
//  Assumes that the first two year characters are '20'.

Function Yymmddhhmmss2time(timeStr)
	String timeStr

	variable  timeValue=date2secs(2000+str2num(timeStr[0,1]), str2num(timeStr[2,3]), str2num(timeStr[4,5])  )
	timeValue+=(3600*str2num(timeStr[6,7])+60*str2num(timeStr[8,9])+str2num(timeStr[10,11])  )
	
	return timeValue

End


// *************************************************
//  Assumes that the first two year characters are 19

Function Yymmddhhmmss2time19(timeStr)
	String timeStr

	variable  timeValue=date2secs(1900+str2num(timeStr[0,1]), str2num(timeStr[2,3]), str2num(timeStr[4,5])  )
	timeValue+=(3600*str2num(timeStr[6,7])+60*str2num(timeStr[8,9])+str2num(timeStr[10,11])  )
	
	return timeValue

End


// *************************************************
// Returns the Igor time value indicated by the string

Function Yyyymmddhhmmss2time(timeStr)
	String timeStr

	variable  timeValue=date2secs(str2num(timeStr[0,3]), str2num(timeStr[4,5]), str2num(timeStr[6,7])  )
	timeValue+=(3600*str2num(timeStr[8,9])+60*str2num(timeStr[10,11])+str2num(timeStr[12,13])  )
	
	return timeValue

End

// New version of ExpandIt handles multiple values by averaging. Also checks for NaNs.
// (Old version only took first data point with multiple time stamps that are equal.)
function ExpandIt(t1, d1, t2, d2, res)
	wave t1, d1, t2, d2
	variable res
	
	variable x1, x2, n1, n2, sumData, nData, startTime, stopTime, firstPoint, firstTime
	variable inc = 0.5*res

	n1 = numpnts(t1)
	n2 = numpnts(t2)
	
	d2 = NaN		//wipe out any values currently in destination wave
	firstPoint = FindFirst(d1, 0)		//find row of first data point
//	if (firstPoint >= n1)
//		Abort "No data found in "+NameOfWave(d1)+" wave. Quiting."
//	endif
	firstTime = t1[firstPoint]		//get time of first data point
	x1 = firstPoint				//index for looping through t1 wave
	startTime = t2[x2] - inc	//0.5*res
	startTime=round(startTime*1000)/1000	//takes care of errors from floating point comparisons of 2 large numbers
	stopTime = t2[x2] + inc	//0.5*res
	stopTime=round(stopTime*1000)/1000
	if (res==0)
		Abort "Resolution can't be 0"
	endif

	if (t2[0] < startTime)			//find row in t2 wave that firstTime is on.
		do
			x2+=1
		while((t2[x2] < startTime) * (x2<n2))
	endif

	do	//loop thru t1 wave
		sumData = 0		//init
		nData = 0		//init

		if ((t1[x1] >= startTime) * (t1[x1] < stopTime))	// & (x1<n1) & (x2<n2)
			do	//loop over res seconds of data
				if (numtype(d1[x1]) == 0)	//only avg non-NaN numbers
					sumData+= d1[x1]	//sum good values
					nData+=1			// number of values in sum
				endif
				x1+=1	//increment row counter for t1 wave
			while((t1[x1] >= startTime) * (t1[x1] < stopTime)  * (x1<n1) * (x2<n2))

			if (nData>0)	//if there was at least 1 value in previous loop
				d2[x2] = (sumData/nData)	//calc avg
//				if (nData>1)	//if there was more than 1 number, print message to user
//					print "Averaging "+num2str(nData)+" points at "+num2str(x1 - 1)
//				endif
			endif
			x2+=1	//increment row counter for t2 wave
//			x1 already incremented above in loop
		else
			if (t1[x1] < startTime)
				x1+=1	//increment index for time1 wave
			else
				x2+=1	//increment index for time2 wave
			endif
		endif
		
		startTime = t2[x2] - inc	//0.5*res
		startTime=round(startTime*1000)/1000	//takes care of errors from floating point comparisons of 2 large numbers
		stopTime = t2[x2] + inc	//0.5*res
		stopTime=round(stopTime*1000)/1000
	while((x1<n1) * (x2<n2))

	if (x2==0)
		Print "*** No times found in common between time waves, so no data was 'expanded' onto new timewave. ***"
	endif
	
end function

// END OF NEW VERSION OF ExpandIt

// *************************************************
// Does the work for many timeexpand macros.
//Note that time waves are truncated to the nearest second.

Function ExpandItOLD(t1,w1,t2,w2)
	Wave t1,w1,t2,w2

	variable x1=0,x2=0
	variable n2=numpnts(t2),n1=numpnts(t1)
	
	t1[,]=round(t1[p])			//*makes timewaves into integers by rounding
	t2[,]=round(t2[p])			//*used to be floor
	
	w2=nan
	
	do
		if ( t2[x2]==t1[x1] )
			w2[x2]=w1[x1]
			x1+=1
			x2+=1
		else
			if ( t2[x2] < t1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while ( ( x2 < n2 )&&(x1< n1 ) )

End


// *************************************************
// Same as ExpandIt except it doesn't round to the second the original time waves

Function ExpandItDontRound(t1,w1,t2,w2)
	Wave t1,w1,t2,w2

	variable x1=0,x2=0
	variable n2=numpnts(t2),n1=numpnts(t1)
	
//	t1[,]=floor(t1[p])			//*makes timewaves into integers by dropping off any fractional seconds
//	t2[,]=floor(t2[p])			//*this is necessary when using Matt W's time values
	
	w2=nan
	
	do
		if ( t2[x2]==t1[x1] )
			w2[x2]=w1[x1]
			x1+=1
			x2+=1
		else
			if ( t2[x2] < t1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while ( ( x2 < n2 )&&(x1< n1 ) )

End


// *************************************************
// Same as ExpandIt but additionally prints the number of non-nans that didn't get 'transferred over' to the newly expanded wave

Function ExpandItPrintMissing(t1,w1,t2,w2, res)
	Wave t1,w1,t2,w2
	variable res

	variable numOrig, numResult
	
 	ExpandIt(t1,w1,t2,w2, res)
 	
	WaveStats/q w1
	numOrig = V_npnts
	
	WaveStats/q w2
	numResult = V_npnts
	
	Print "// There are "+num2str(numOrig - numResult) +" number of values out of "+num2str(numOrig)+" not represented in "+NameofWave(w2) +" after expanding it. "

End


//*************************************
// Expands data wave when original data wave may have duplicate times and occasional nans.
// For example, if original time, data pair has consequtive values such as
// 12:00:00 (nan)
// 12:00:00 4
// This expandIt will not transfer over the nan, it will use the first non-nan, in this case, 4.

Function ExpandItRoundTimePossibleDup(t1,w1,t2,w2)
	Wave t1,w1,t2,w2

	variable x1=0,x2=0
	variable n2=numpnts(t2),n1=numpnts(t1)
	
	t1[,]=round(t1[p])	
	t2[,]=round(t2[p])			
	
	w2=nan
	
	do
		if ( t2[x2]==t1[x1] )
			if ( numtype(w1[x1])==0 )
				w2[x2]=w1[x1]
				x1+=1
				x2+=1
			else
				x1+=1
			endif
			
		else
			if ( t2[x2] < t1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while ( ( x2 < n2 )&&(x1< n1 ) )
end

// *************************************************
// Performs an expansion similar to the ExpandIt function, only is used for faster than 1 hz data.

Function ExpandItxHz(t1,w1,t2,w2, xHzval)
	Wave t1,w1,t2,w2
	Variable xHzval

	variable x1=0,x2=0
	variable n2=numpnts(t2),n1=numpnts(t1)
	
	t1[,]=floor(t1[p]*xHzval)/xHzval			//* make sure time is in even fractional seconds
	t2[,]=floor(t2[p]*xHzval)/xHzval			
	
	w2=nan
	
	do
		if ( t2[x2]==t1[x1] )
			w2[x2]=w1[x1]
			x1+=1
			x2+=1
		else
			if ( t2[x2] < t1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while ( ( x2 < n2 )&&(x1< n1 ) )

End

// *************************************************
// Performs an expansion but doesn't round t1 and t2 to even fractional seconds. Finds nearest time in t2 within the supplied delta.
// Especially useful for data at faster than 1Hz that isn't on an even time base.

function ExpandItFindNearest(t1, w1, t2, w2, delta)
	wave t1, w1, t2, w2
	variable delta
	
	variable x1, x2, n1, n2 , diff1, diff2
	n1 = numpnts(t1)	//length of timewave1
	n2 = numpnts(t2)	//length of timewave2

	do
		diff1 = abs(t2[x2] - t1[x1])
		diff2 = abs(t2[x2] - t1[x1+1])
		if ((diff1 <=  diff2) * (diff1 <= delta) * (numtype(w1[x1]) == 0))    //if this difference is smaller than the next difference, and difference is less than delta, and isn't NaN.
			w2[x2]=w1[x1]
			x1+=1
			if ((abs(t1[x1] - t2[x2])) < (abs(t1[x1+1] - t2[x2])))	//only increment x2 if the NEXT point in t1 is farther from t2 than the current point in t1
				x2+=1
			endif
		else
			if ( t2[x2] < t1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while((x2 < n2 )&&(x1< n1))
end function

// *************************************************
// Expands data which is on StartTime/StopTime, onto another wave such as AOCTimewave, and fills in the values between
//	the start/stop time with the value that was in the data wave.
// Useful for example to expand CH2O data from Start/Stop time onto AOCtimewave. 
Function ExpandItFill(start1, stop1, w1, t2, w2)
	Wave start1, stop1, w1, t2, w2

	variable x1=0, x2=0
	variable n2=numpnts(t2), n1=numpnts(start1)
		
	w2=nan
	
	do
		if (( t2[x2] >= start1[x1] ) * ( t2[x2] < stop1[x1] ))	//need to fill in between values
			do
				w2[x2] = w1[x1]	//set to same value between start and stoptime
				x2+=1
				if (x2>n2)	//exceeded length of wave t2, so there is no more time to expand onto or fill
					break
				endif
			while(( t2[x2] >= start1[x1] ) * ( t2[x2] < stop1[x1] ))
			w2[x2] = w1[x1]
			x1+=1
			x2+=1
		else
			if ( t2[x2] < start1[x1] )
				x2+=1
			else
				x1+=1
			endif
		endif
	while ( ( x2 < n2 )&&(x1< n1 ) )

End

// *************************************************
// "Expands" data in wave1 by duplicating points where time waves overlap and puts nans elsewhere. 
// 	A new wave is created from the user input string (given at the prompt).
//	If wave1 is double precision, wave2 will be double also. Else, the default is single.

Macro TimeExpandCreateWave(time1, wave1, time2, wave2, res)
	String time1, wave1, time2, wave2="new_expanded"
	variable res=1 //# of seconds to average to
	prompt time1,"Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Name of resulting expanded data wave?"
	prompt res, "Enter number of seconds to average:"
	
	Silent 1
	
	if (numpnts($time1)!=numpnts($wave1) ) )
	Print "// ", time1, wave1, numpnts($time1), numpnts($wave1)
	Abort "The number of points must be the same; see the history window to help debug. Aborting from TimeExpandCreateWave"
	endif
	
	wave2 = CleanUpName(wave2,0)
	Variable precision
	precision = wavetype($wave1)
	if (precision==4)
	Make/O/D/N=(numpnts($time2)) $wave2 // this is a double precision wave
	else
	Make/o/d/n=(numpnts($time2)) $wave2 // this is a single precision wave
	endif
	$wave2=nan
	
	ExpandIt($time1,$wave1,$time2,$wave2, res)

End

// *************************************************
// "Expands" data values in wave1 by duplicating points where time waves overlap and puts nans elsewhere.
// The only wave which is altered in this macro is wave2.

Macro TimeExpandDontRoundCreateWave(time1, wave1, time2, wave2)
	String time1, wave1, time2, wave2="new_expanded"
	prompt time1,"Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Name of resulting expanded data wave?"

	Silent 1
	
	if (  numpnts($time1)!=numpnts($wave1) ) )
		Print "//", time1, wave1, numpnts($time1), numpnts($wave1) 
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from TimeExpandDontRoundCreateWave"
	endif
	
	wave2 = CleanUpName(wave2,0)
	Make/o/d/n=(numpnts($time2)) $wave2
	$wave2=nan
	
	ExpandItDontRound($time1,$wave1,$time2,$wave2, res)

End


// *************************************************
// "Expands" data in wave1 by duplicating points where time waves overlap and puts nans elsewhere
//	The only wave which is altered in this macro is wave2

Macro TimeExpandDontRoundTime(time1, wave1, time2, wave2)
	String time1, wave1, time2, wave2="new_expanded"
	prompt time1,"Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)

	Silent 1
	
	if (  numpnts($time1)!=numpnts($wave1) ) )
		Print "//", time1, wave1, numpnts($time1), numpnts($wave1) 
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from TimeExpandDontRoundTime"
	endif
	
	wave2 = CleanUpName(wave2,0)
	Make/o/d/n=(numpnts($time2)) $wave2
	$wave2=nan
	
	ExpandItDontRound($time1,$wave1,$time2,$wave2)

End


// *************************************************
// Same as Time_expand except it prints the number of points that did not get "copied over" to the expanded wave.
// This is useful for determining if the original timewaves had duplicate values for one second
// that did not get transferred onto the expanded wave.

Macro TimeExpandPrintMissCreateWave(time1, wave1, time2, wave2, res)
	String time1, wave1, time2, wave2
	variable res=1
	prompt time1, "Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt res, "Number of seconds to average:"

	Silent 1
	
	if (  numpnts($time1)!=numpnts($wave1) )||( numpnts($time2)!=numpnts($wave2) )
		if (numpnts($time1)!=numpnts($wave1) )
			Print "//", time1, wave1, numpnts($time1), numpnts($wave1) 
		else
			Print "//", time2 wave2, numpnts($time2), numpnts($wave2) 
		endif
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from TimeExpandPrintMissCreateWave"
	endif
	
	$wave2=nan
	
	ExpandItPrintMissing($time1,$wave1,$time2,$wave2, res)

End


// *************************************************
// Same as Time_expand except it prints the number of points that did not get  "copied over" to the expanded wave.
// This is useful for determining if the original timewaves had duplicate values for one second
// that did not get represented, or transfered over to the newly expanded wave.

Macro TimeExpandPrintMissing(time1, wave1, time2, wave2, res)
	String time1, wave1, time2, wave2
	variable res=1
	prompt time1, "Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt res, "Number of seconds to average:"

	silent 1
	
	if (numpnts($time1)!=numpnts($wave1) )
		Print "// ", time1, wave1, numpnts($time1), numpnts($wave1) 
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from TimeExpandPrintMissing"
	endif
	
	wave2 = CleanUpName(wave2,0)
	Make/o/d/n=(numpnts($time2)) $wave2
	$wave2=nan
		
	ExpandItPrintMissing($time1,$wave1,$time2,$wave2)

End


// *************************************************
// "Expands" data in wave1 by duplicating points where time waves overlap and puts nans elsewhere. 

Macro Time_expand(time1, wave1, time2, wave2, res)
	String time1, wave1, time2, wave2
	variable res=1
	prompt time1, "Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""),";",4)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt res, "Number of seconds to average:"

	silent 1
	
	if (numpnts($time1)!=numpnts($wave1) )||( numpnts($time2)!=numpnts($wave2))
		if (numpnts($time1)!=numpnts($wave1))
			Print "// ", time1, wave1, numpnts($time1), numpnts($wave1) 
		else
			Print "// ", time2 wave2, numpnts($time2), numpnts($wave2) 
		endif
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from Time_expand"
	endif
	
	if (res==0)
		Abort "Resolution can't be zero."
	endif
	
	$wave2=nan
	
	ExpandIt($time1,$wave1,$time2,$wave2, res)

End


// *************************************************
// "Expands" data in wave1 by finding closest  point in time2 wave within the specified delta. 
//	No waves are rounded off. The destination wave, wave2 is NaNed before this function executes.
// Useful for expanding 1 Hz data onto faster data that isn't on regular interval and not on even increments.

Macro TimeExpandFindNearest(time1, wave1, time2, wave2, delta)
	String time1, wave1, time2, wave2
	variable delta
	prompt time1, "Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt delta, "Maximum delta from nearest value?"

	silent 1
	
	if (numpnts($time1)!=numpnts($wave1) )||( numpnts($time2)!=numpnts($wave2))
		if (numpnts($time1)!=numpnts($wave1))
			Print "// ", time1, wave1, numpnts($time1), numpnts($wave1) 
		else
			Print "// ", time2 wave2, numpnts($time2), numpnts($wave2) 
		endif
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from Time_expand"
	endif
	
	$wave2=nan
	
	ExpandItFindNearest($time1, $wave1, $time2, $wave2, delta)

End

// *************************************************
// "Expands" data in wave1 by duplicating points where time waves overlap and puts nans elsewhere. 
//	Note that time waves are truncated to the nearest second (in the ExpandIt function)

Macro TimeExpandFill(starttime1, stoptime1, wave1, time2, wave2)
	String starttime1, stoptime1, wave1, time2, wave2
	prompt starttime1, "Start Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt stoptime1, "Stop Time wave for the data you want to expand?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave1, "Data wave you want expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt time2, "Timewave to which data will be expanded?",popup SortList(WaveList("*",";",""), ";", 16)
	prompt wave2, "Resulting expanded data wave?",popup SortList(WaveList("*",";",""), ";", 16)

	silent 1
	
	if (numpnts($starttime1)!=numpnts($wave1) )|| (numpnts($stoptime1)!=numpnts($wave1) ) || ( numpnts($time2)!=numpnts($wave2))
		if (numpnts($starttime1)!=numpnts($wave1))
			Print "// ", starttime1, wave1, numpnts($starttime1), numpnts($wave1) 
		else
			Print "// ", time2 wave2, numpnts($time2), numpnts($wave2) 
		endif
		Abort "The number of points must be the same; see the history window to help debug.  - Aborting from Time_Expand_Fill"
	endif
	
	$wave2=nan
	
	ExpandItFill($starttime1, $stoptime1, $wave1, $time2, $wave2)

End


// *************************************************
// Kills everything but notebooks, procedure files and xop windows.

Function KillAbsolutelyEverything()

	KillAllTables()
	KillAllGraphs()
	KillAllLayouts()
	KillAllPanels()
	KillPath/A/Z
	
	SetDataFolder root:
	KillDataFolder root:		//kills waves, strings, variables, child data folders.

End


// *************************************************
//  Kills all graphs (duh).

Function KillAllGraphs()

	KillTheseWindows(1)

End


// *************************************************
//  Kills all layouts (duh).

Function KillAllLayouts()

	KillTheseWindows(4)

End


// *************************************************
//  Kills all panels (duh).

Function KillAllPanels()

	KillTheseWindows(64)

End


// *************************************************
//  Kills all tables (duh).

Function KillAllTables()

	KillTheseWindows(2)

End


// *************************************************
//  Kills all tables graphs layouts and waves (duh).

Function KillAllTablesGraphsLayoutsWaves()

	KillAllTables()
	KillAllGraphs()
	KillAllLayouts()
	KillWaves/A/Z		// kills waves only in current data folder

End


// *************************************************
// Kills the windows specified by the igor code winTypeNum

Function KillTheseWindows(winTypeNum)
	Variable winTypeNum

	string cmd 
	
	if ( strlen(WinName(0, winTypeNum))>0 )	
		do
			cmd =  "DoWindow/K "+WinName(0,winTypeNum)
			Execute cmd			
		while(strlen(WinName(0,winTypeNum))>0)		
	endif

End


// *************************************************
//  Kills all killable waves in list.
//  It does not generate an error if one of the waves in the list is in use (and cannot be deleted)
//  It does not generate an error if the wave in the list does not exist
// It also used to be called KillTheseWaves

Function KillWavesInList(wList)
	String wList

	variable oneWaveVar = 0, oneWaveNum = ItemsInList(wList)
	string oneWaveStr 
	
	for (oneWaveVar=0; oneWaveVar<oneWaveNum; oneWaveVar+=1)
		oneWaveStr = StringFromList(oneWaveVar, wList, ";")
		if (exists(oneWaveStr)==1)
			Killwaves/Z $oneWaveStr
		else
	//		print  "This one doesn't exist "+oneWaveStr
		endif
	endfor

End


// *************************************************
//  Kills all killable zero-length waves in list.
//  It does not generate an error if one of the waves in the list is in use (and cannot be deleted)
//  It does not generate an error if the wave in the list does not exist

Function KillZeroLengthWaves(wList)
	String wList

	variable oneWaveVar = 0, oneWaveNum = ItemsInList(wList)
	string oneWaveStr 
	
	for (oneWaveVar=0; 	oneWaveVar<oneWaveNum; oneWaveVar+=1)
		oneWaveStr = StringFromList(oneWaveVar, wList, ";")
		if (exists(oneWaveStr)==1)
			if (numpnts($oneWaveStr)==0)
				Killwaves/Z $oneWaveStr
			endif
		else
			print  "// This one doesn't exist "+oneWaveStr+"; printing from KillZeroLengthWaves"
		endif
	endfor
	
End


// *************************************************
// Checks to see if all the waves with name entries in a list are of the same length.
// Returns true, 1, if all waves in list are same length or if string is empty.

Function AllInListSameLength(wList)
	String wList

	string wStr
	variable lvar, numInWave, returnVal = 1

	wStr = StringFromList(lvar, wList, ";")
	numInWave = numpnts($(wStr))
	
	for(lvar=0; lvar<ItemsInList(wList); lvar+=1)
		wStr =  StringFromList(lvar, wList, ";")
		if (Exists(wStr)==1 )
			if (numInWave != numpnts($(wStr))  )
				returnVal = 0
			endif
		endif
	endfor

	return returnVal
	
End


// *************************************************
//Returns a string from the list that looks like str1 but not str2

Function/s FindStr1NotStr2InList(list, str1, str2)
	String list, str1, str2

	string str, returnStr = ""
	variable i, n = ItemsInList(list)
	
	for (i = 0; i <n; i+=1)
		str = StringFromList(i, list, ";")
		if (  (strsearch(str, str1, 0 )>0)&& (strsearch(str, str2, 0 )==-1 )  )
			returnStr = str
		endif
	endfor
	
	return returnStr

End


// *************************************************
//  Returns the name of a text wave containing elements in the list; overwrites contents of text wave textWaveNameStr if it exists

Function/S List2textWave(ListStr, textWaveNameStr)
	String ListStr, textWaveNameStr

	string returnStr = textWaveNameStr
	variable i=0, n=ItemsinList(ListStr)
	
//	if (strlen(ListStr)==0)
//		abort "There were no items in this list"
//	endif
	
	if (Exists(textWaveNameStr)==1 )
		WAVE/T tw = $textWaveNameStr
		Redimension/N=(n)  tw
	else
		returnStr = CleanupName(textWaveNameStr, 0 )
		if (CheckName(returnStr, 1)!=0)
			 returnStr = UniqueName(textWaveNameStr[0, min(20, strlen(returnStr)-1)], 1, 0 )
		endif
		Make/o/t/n=(n) $returnStr
		WAVE/T tw = $returnStr
	endif
	
	for (i=0; i<n; i+=1)
		tw[i]  = StringFromList(i, ListStr, ";")
	endfor
	
	return returnStr
	
End


// *************************************************
// Same as List2textWave except the function is supplied with the text wave, not the name of the text wave to be created.

Function List2textWaveNoCreate(listStr, textWaveName)
	String listStr
	Wave/t textWaveName

	variable i=0, n=ItemsinList(ListStr)
	
//	if (strlen(ListStr)==0)
//		abort "There were no items in this list"
//	endif
			
	if (numpnts(textWaveName) != n)
		Abort "The number of points in the wave "+num2str(numpnts(textWaveName))+" is not the same as the list "+num2str(n)+"List2textWaveNoCreate - Aborting from List2textWaveNoCreate"
	endif
	
	for (i=0; i<n; i+=1)
		textWaveName[i]  = StringFromList(i, listStr, ";")
	endfor

End


// *************************************************
// Joins 2 lists without duplication of entries of the second list into the first.

Function/s MergeLists(list1, list2)
	String list1, list2

	string str2, returnStr = list1
	
	variable i, n= ItemsInList(list2)
	for(i=0 ; i < n ; i+=1)	
		str2 = StringFromList(i, list2,";")
		if (WhichListItem(str2 , returnStr, ";" )==-1)
			returnStr+=str2+";"
		endif					
	endfor										
	
	return returnStr
	
End


// *************************************************
// Returns a list with every item in removeList removed.

Function/s RemoveItemsInListFromList(removeList, masterList)
	String removelist, masterList

	string newMasterList = masterList
	variable i, n
	
	n=itemsinlist(removelist)
	
	for(i=0;i<n;i+=1)
		newMasterList = RemoveFromList(StringFromList(i, removeList, ";"), newMasterList, ";")
	endfor
	
	return newMasterList

End


// *************************************************
// Returns a string containing the common entries in the two lists.  (Sort of a list AND)

Function/s SameInBothLists(list1, list2)
	String list1, list2

	string str2, returnStr =""
	
	variable i, n= ItemsInList(list1)
	for(i=0; i < n; i+=1)	
		str2 = StringFromList(i, list1,";")
		if (WhichListItem(str2 , list2, ";" )>=0)
			returnStr+=str2+";"
		endif					
	endfor										
	
	return returnStr
	
End


// *************************************************
// Sorts a list alphabetically - contribution from Harald

Function/S SortListAlphabetically(listStr)
	String listStr

	string textWav, returnStr
	
	textWav=UniqueName("textWave",1,0)
	Make/O/T/n=(ItemsInList(listStr)) $textWav
	
	textWav = List2textWave(listStr, textWav)
	Sort $textWav $textWav //sort textWave alphabetically
	
	returnStr = textWave2List($textWav)
	Killwaves $textWav
	
	return returnStr

End


// *************************************************
// Returns a string containing a ; separated list of text wave contents.

Function/s TextWave2List(textWave)
	Wave/t textWave

	string returnList = ""
	
	variable i=0, n = numpnts(textWave)
	
	if ( n==0 )
		Print "// There were no items in this textWave "+ NameOfWave(textWave) 
	endif
	
	for (i=0; i<n; i+=1)
		returnList+=textWave[i]+";"
	endfor
	
	return returnList

End


// *************************************************
// Returns a string list containing the names of all text waves with names matching matchstr.

Function/s TextWaveList(matchStr,separatorStr,optionsStr)
	String matchStr, separatorStr, optionsStr

	string waveListStart = WaveList(matchStr,separatorStr,optionsStr)
	string waveListReturn = ""
	
	variable waveListIndex=0, numItems = ItemsInList(waveListStart)
	string waveListStr 
	
	for(waveListIndex=0;waveListIndex<numItems; waveListIndex+=1)
		waveListStr = StringFromList(waveListIndex, waveListStart, ";")
		if(NumberByKey( "NUMTYPE", waveinfo($waveListStr,0) )==0)
			waveListReturn+=waveListStr+";"
		endif
	endfor
	
	return waveListReturn

End


// *************************************************
// Returns a formatted string of the OS path of PathNameStr whereby
//  each folder is on a separate line and is indented one tab from the parent directory.
//  Example: suppose a path is created as follows....
//  NewPath testpath "Macintosh HD:Users:sueper:Documents:datahere:"
//  Print convertpathStr("testpath") results in 
//  Macintosh HD
//	Users
//		sueper
//			Documents
//				datahere
// This string can then be nicely uesd for titleboxes 
 
Function/s ConvertPathStr(pathNameStr)
	String pathNameStr

	string returnStr
	
	string folderStr, tabStr
	variable i

	returnStr = ""
	tabStr = ""
	PathInfo $pathNameStr
	
	if (strlen(S_path)>0)	
		do
			folderStr = ParseFilePath(0, S_path, ":", 0, i)
			if (strlen(folderStr))
				returnStr+=tabStr+folderStr+"\r"
			endif
			tabStr+="\t"
				
			i+=1
		while(strlen(folderStr)>0)
		returnStr = returnStr[0,strlen(returnStr)-2]
	endif

	return returnStr

End


// *************************************************
// Creates a new directory and igor path to this directory which will be on the same level (a sister folder) 
// to the directory where OldPathStr points to.

Function CreateANewFolderPathFromPath(oldPathStr,  newFolderStr, newPathStr)
	String oldPathStr,  newFolderStr, newPathStr

	variable selectNewFoldersVar = 1

	Pathinfo/S $oldPathStr

	if (V_Flag ==0)
		Abort "There is no path with the name "+oldPathStr+" - Aborting from CreateaNewFolderPathFromPath"
	endif
	
	string rootPathStr = GetPathRootString(S_path)
	string dummyStr = rootPathStr
	
	prompt selectNewFoldersVar, "Enter 0 to escape, Enter 1 to create a folder with the following name at the location below"
	prompt dummyStr, "The location where the folder that will be created.: (DO NOT EDIT)"
	prompt newFolderStr, "Enter the name of the new folder to be created and placed in the folder identified above."
	
	doPrompt "Enter the values for creating a new folder.", selectNewFoldersVar, dummyStr, newFolderStr

	if (selectNewFoldersVar == 1)	// we did not escape
		NewPath/C/O $newPathStr rootPathStr+CleanupName(newFolderStr, 0)
	endif
	
End


// *************************************************
// Returns a string containing the name of the directory where pathStr resides.
// This is intended to be used as a means of separating an upperlevel directory from a path string
// It is assumed that the parameter passed is a string such as "Nice HD #2:convert ds files:"
// which can be obtained from a Pathinfo call
//	The code has been replaced with a call to ParseFilePath().

// ParentFolder returns "HD:folder:" from something like "HD:folder:subfolder:theFile.ext"

Function/s GetPathRootString(pathStr)
	String pathStr

	return ParseFilePath(1, pathStr, ":", 1, 1)

End


// *************************************************
// Loads waves from an igor path.
// If you don't send this function a string,  a dialog will prompt for a path
// All current igor waves will be overwritten.

Function LoadIgorWavesFromThisPath(pathStr)
	String pathStr
	Prompt pathStr,"From which path would you like to load igor waves?",popup PathList("*",";","")

	PathInfo $PathStr
	
	if(V_flag==0)
		Abort "There is no such path - Aborting from LoadIgorWavesFromThisPath"
	endif
	
	variable fileIndex
	string fileName
	
	do
		fileName= IndexedFile($PathStr,fileIndex,".ibw")
		if (strlen(fileName) != 0)
	      		LoadWave/o/H/Q/P=$PathStr fileName
	      endif
	      fileIndex += 1
	while (strlen(fileName) !=0)
	
	Print "// ", num2str(fileIndex-1) +" many waves have been loaded."

End


// *************************************************
// Puts the name of the files in the text wave fileNameTextWave which exist in the path defined by pathStr.

Function PopulateFileNameTextWave(fileNameTextWave, pathStr)
	Wave/t fileNametextWave
	String pathStr

	string wlist, thisfile
	variable i, num
	
	Redimension/N=0 fileNametextWave

	PathInfo $pathStr
	
	if (V_Flag == 0)
		NewPath/q/o /m="Where is the folder with the files?"  $pathStr
	endif

	wlist = IndexedFile($pathStr, -1, "????" )
	num = ItemsInList(wlist, ";")
	for(i=0; i<num; i+=1)
		thisfile = StringFromList(i, wlist, ";")
		insertPoints i,  1, fileNameTextWave
		fileNameTextWave[i] = thisfile
	endfor
End

// *************************************************
// Checks to see if a path exists; if it does, returns 1
//  If the path does not exist, it creates it, prompts the user for the path. returns 1 if successful, 0 if not.

Function SetUpPath(pathNameStr, messageStr)
	String pathNameStr, messageStr

	variable returnVal = 1
	
	PathInfo $pathNameStr
	
	if(V_flag==0)
		NewPath/m=messageStr $pathNameStr
		if (V_flag==0)
			returnVal = 1
		else 
			returnVal = 0
		endif
	else
		returnVal = 1
	endif

	return returnVal

End


// *************************************************
// Adds a suffix to all wave names in the list

Function AppendSuffixToWaveList(replaceList, suffix)
	String replaceList, suffix

	string replaceOne
	variable i=0, n = ItemsInList(replaceList)

	for(i=0;i<n; i+=1)
		replaceOne = StringFromList(i, replaceList, ";")
		if (CheckName( (replaceOne[0,strlen(replaceOne)-1]+suffix), 1)==0)  // valid name
			Rename $replaceOne $(replaceOne[0,strlen(replaceOne)-1]+suffix)
		else
			Rename $replaceOne $(UniqueName((replaceOne[0,strlen(replaceOne)-1]+suffix), 1, 0)  )
		endif		
	endfor

End


// *************************************************
//Does the same as RenameWavesWithThisSuffix except that it doesn't search all waves just the waves in wList. 
 
Function RenameWaveListWithThisSuffix(wList, suffix2replace, newSuffix)
	String wList, suffix2replace, newSuffix

	string replaceOne, replaceList =  SameInBothLists(wList, Wavelist("*"+suffix2replace, ";", "")  )
	variable strPlace, i=0, n = ItemsInList(replaceList)

	for(i=0;i<n; i+=1)
		replaceOne = StringFromList(i, replaceList, ";")
		strPlace = StrLen(replaceOne)
		if (CheckName(replaceOne[0,strPlace-1]+newSuffix, 1)==0)  // valild name
			Rename $replaceOne $(replaceOne[0,strPlace-1]+newSuffix)
		else
			Rename $replaceOne $(UniqueName(replaceOne[0,strPlace-1]+newSuffix, 1, 0)  )
		endif		
	endfor

End


// *************************************************
// Renames waves beginning with prefix2replace in current data folder to have a new prefix.

Function RenameWavesWithThisPrefix(prefix2replace, newprefix)
	String prefix2replace, newprefix

	string replaceOne, replaceList = Wavelist(prefix2replace+"*", ";", "")
	variable strPlace, i=0, n = ItemsInList(replaceList)

	for(i=0;i<n; i+=1)
		replaceOne = StringFromList(i, replaceList, ";")
		strPlace = StrLen(prefix2replace)
		if (CheckName(newprefix+replaceOne[strPlace,StrLen(replaceOne)-1], 1)==0)  // valild name
			Rename $replaceOne $(newprefix+replaceOne[strPlace,StrLen(replaceOne)-1])
		else
			Rename $replaceOne $(UniqueName(newprefix+replaceOne[strPlace,StrLen(replaceOne)-1], 1, 0)  )
			Print "// printing from RenameWavesWithThisPrefix; Badname ", replaceOne, UniqueName(newprefix+replaceOne[strPlace,StrLen(replaceOne)-1], 1,0)			
		endif		
	endfor

End


// *************************************************
// Searches all waves with suffix2replace.

Function RenameWavesWithThisSuffix(suffix2replace, newSuffix)
	String suffix2replace, newSuffix

	string replaceOne, replaceList = Wavelist("*"+suffix2replace, ";", "")
	variable strPlace, i=0, n = ItemsInList(replaceList)

	for(i=0;i<n; i+=1)
		replaceOne = StringFromList(i, replaceList, ";")
		strPlace = StrSearch(replaceOne, suffix2replace,0)
		if (CheckName(replaceOne[0,strPlace-1]+newSuffix, 1)==0)  // valild name
			Rename $replaceOne $(replaceOne[0,strPlace-1]+newSuffix)
		else
			Rename $replaceOne $(UniqueName(replaceOne[0,strPlace-1]+newSuffix, 1, 0)  )
		endif		
	endfor

End


// *************************************************
// Saves one wave as a text file with format as given in formatStr.

Function SaveWaveAsText(pathStr, waveNameStr, formatStr)
	String pathStr, waveNameStr, formatStr

	Pathinfo $pathStr
	if (V_flag == 0)
		NewPath/m="Where do you want to put the text files?" $pathStr
		Pathinfo $pathStr
	endif
	
	string tempWaveStr, cmd
	variable fileNum
	
	tempWaveStr= UniqueName("temporary_wave", 0, 0)
	
	Duplicate/o $waveNameStr $tempWaveStr
//	ChangeNans($tempWaveStr,-999)

	formatStr+="\\r\\n"
	
	Open/p=$pathStr fileNum as waveNameStr+".txt"
	if (strlen(S_fileName)>0)
		Fprintf fileNum, waveNameStr+"\r\n"
		cmd = "Wfprintf  "+num2str(fileNum)+", \""+ formatStr+"\" "+tempWaveStr 
	 
		Execute cmd
		Close fileNum
	endif
	
	Killwaves $tempWaveStr
	
End

		
// *************************************************
// Saves waves in List as igor binary files

Function SaveWavesInList(pathStr, waveListStr)
	String pathStr, waveListStr
		
	variable oneWaveVar = 0, oneWaveNum = ItemsInList(waveListStr)
	string oneWaveStr
	
	Pathinfo $pathStr
	if (V_flag == 0)
		NewPath/m="Where do you want to put the igor files?" $pathStr
		Pathinfo $pathStr
	endif

	for (oneWaveVar=0; 	oneWaveVar<oneWaveNum; oneWaveVar+=1)
		oneWaveStr = StringFromList(oneWaveVar, waveListStr, ";")
		if (exists(oneWaveStr)==1)
			Save/C/O/P=$(pathStr) $oneWaveStr as oneWaveStr+".ibw"		
		else
			print  "// Couldn't save this one because it doesn't exist: "+oneWaveStr+"; printing from SaveWavesInList"
		endif
	endfor
		
End


// *************************************************
// Saves waves in List as text files
// Originally function SaveWavelistAsText(waveListStr, pathStr,  formattype)

Function SaveWavesInListAsText(pathStr, waveListStr,  formattype)
	String pathStr, waveListStr, formattype
		
	variable oneWaveVar = 0, oneWaveNum = ItemsInList(waveListStr)
	string oneWaveStr
	
	Pathinfo $pathStr
	if (V_flag == 0)
		NewPath/m="Where do you want to put the text files?" $pathStr
		Pathinfo $pathStr
	endif

	for(oneWaveVar=0; oneWaveVar<oneWaveNum; oneWaveVar+=1)
		oneWaveStr = StringFromList( oneWaveVar, waveListStr,";")
		Execute "SaveWaveAsText(\""+pathStr+"\", \""+oneWaveStr+"\", \""+formattype+"\")"
	endfor
		
End


// *************************************************
// Transforms the textwave to a string then calls SaveWavesInList

Function SaveWavesInTextWave(pathStr, textWave)
	String pathStr
	Wave/t textWave	

	SaveWavesInList(pathStr, textWave2List(textWave))
		
End


// *************************************************
// Returns the number of carriage returns in a string

Function CountCarriageReturns(str)
	String str

	variable i = 0, n = 0
	
	do
		i = strsearch(str, "\r", i)
		if (i>=0)
			n+=1
			i+=1
		endif
	while(i>=0)
	
	return n

End


// *************************************************
// Returns a string which has the same point number as the
// first entry in valWave which has the same value as val2match
//  If no match was found, the function returns "NoMatchFound"
// If more than one match, it gives the last value.

Function/s FindValMatchText(val2match, valWave, textWave)
	Variable val2match
	Wave valWave
	Wave/t textWave

	variable i=0, n=numpnts(valWave)
	string returnStr="NoMatchFound"
	
	if ( numpnts(textWave) != n )
		Abort "There must be the same number of points in both waves - Aborting from FindValMatchText"
	endif
	 
	do
		if (val2match==valWave[i])
	 		returnStr = textWave[i]
	 		i=n
	 	endif
	 
	 	i+=1
	while(i<n) 
	 
	return returnStr
 
End
 

// *************************************************
// Appends the padStr to items in the textwave

Function FrontPadTextWave(textWave, padStr, padLength)
	Wave/t textWave
	String padStr
	Variable padLength

	variable i=0, n=numpnts(textWave)
	
	do
		if (strlen(textWave[i]) <padLength)
			textWave[i] = padStr + textWave[i]
		endif
		i+=1
	while(i<n)

End


// *************************************************

// Prints a string to a file.  If the string doesn't have  line feeds or carriage returns after 1000 chars
// this function will insert some so that we don't get any errors. 
// Donna added a noCRflag optional flag for when you want the string to just print.
Function FPrintFLongerThan1000Chars(fileNum, possiblyLongStr, [noCRflag])
	Variable fileNum, noCRflag
	String PossiblyLongStr

	string str, CRstr = "\r\n"
	variable pos1, pos2, n = strlen(possiblyLongStr)
	
	pos1 = 0
	pos2 = min(pos1+999, n-1)
	
	if (! ParamIsDefault(noCRflag) && noCRflag==1)
		CRstr=""
	endif
	
	do
		str = possiblyLongStr[pos1, pos2]
		if(strlen(str)<1000)	
			Fprintf fileNum, "%s"+CRstr, str
		elseif ( (strsearch(str, "\r", 0)<0) && (strsearch(str, "\n", 0)<0) )   	// force inserting some end of CRLFs
			Fprintf fileNum, "%s"+CRstr, str
		else
			Fprintf fileNum, "%s", str
		endif
		pos1 = pos2+1
		pos2 = min(pos1+999, n-1)
	while(pos1<n-1)
	
End


// *************************************************
//  Returns a 2 char or more length string of the number x

Function/s Num2Str2Char(i)
	Variable i

	string returnStr = num2str(i)
	
	if ( (i>=0) && (i<=9) )
		returnStr = "0"+returnStr
	endif
	
	return returnStr

End


// *************************************************
//  Returns a 3 char or more length string of the number x

Function/s Num2Str3Char(i)
	Variable i

	string returnStr = num2str(i)
	
	if ( (i>=0) && (i<=9) )
		returnStr = "00"+returnStr
	else
		if( (i>=10) && (i<=99) )
			returnStr = "0"+returnStr		
		endif	
	endif
	
	return returnStr

End


// *************************************************
// Returns a string whereby the initial carriage returns are removed but remaining line feeds are replaced with one space.

Function/s RemoveCarriageReturns(str)
	String str
	
	return ReplaceTheseCharsWithSpace("\r", str)

End


// *************************************************
// Returns a string with all double spaces in str removed; is iterative.

Function/s RemoveDoubleSpaces(str)
	String str
	
	string returnStr = str
	variable n
	
	do
		n = strlen (returnStr)
		returnStr =ReplaceString("  ",returnStr," ")   //  StrSubstitute("  ",returnStr," ")
	while(n!=strlen(returnStr) )
	
	return returnStr

End


// *************************************************
//  Redimensions tw so that it doesn't have an empty string in any rows.

Function RemoveEmptyStringFromTextWave(tw)
	Wave/t tw

	variable p, numPoints, numEmptyStrs
	string str
	
	numEmptyStrs = 0
	p = 0										// the loop index
	numPoints = numpnts(tw)						// number of times to loop

	do
		str = tw[p]
		if (strlen(str)==0)					
			numEmptyStrs += 1
		else										
			tw[p - numEmptyStrs] = str			// copy to input wave
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numEmptyStrs, numEmptyStrs, tw
	
	return(numEmptyStrs)
End


// *************************************************
// Returns a string with all  preceding spaces Iteratively removed.

Function/s RemoveFirstSpaces(str)
	String str

	variable i=0, n, isSpace
	string returnStr = str
	
	do
		n = strlen(returnStr)
		if (n!=0)
			isSpace = cmpstr(returnStr[0], " ")
			if (isSpace==0)
				returnStr = returnStr[1, n-1]
			endif
		endif
	while(  (isSpace==0) &&(n!=0)  )
	
	return returnStr
End


// *************************************************
// Returns a string whereby the Initial line feeds are removed but remaining line feeds are replaced with one space.

Function/s RemoveLineFeeds(str)
	String str
	
	return ReplaceTheseCharsWithSpace("\n", str)

End


// *************************************************
// Returns a string whereby the initial tabs are removed but remaining line feeds are replaced with one space.

Function/s RemoveTabs(str)
	String str

	return ReplaceTheseCharsWithSpace("\t", str)

End


// *************************************************
// Same as RemoveFirstSpaces, except for any char string instead of spaces.

Function/s RemoveTheseFirstChar(str, charStr)
	String str, charStr

	variable i=0, n, isSame
	string returnStr = str
	string thischar

	n = strlen(str)
	
	if (n==0)
		return ""
	endif
	
	do
		isSame = cmpstr(str[i], charStr[i])
		if (isSame==0)
			returnStr = str[i+1, n-1]
		endif
		i+=1
	while(  (isSame==0) &&(n!=0)  )
	
	return returnStr

End


// *************************************************
// Returns a string whereby carriage returns are removed from the end of the string

Function/s RemoveTrailingCarriageReturns(str)
	String str
	
	return RemoveTrailingChars("\r", str)

End


// *************************************************
// Returns a string whereby carriage return -  line feeds are removed from the end of the string

Function/s RemoveTrailingChars(char2remove, str)
	String char2remove
	String str

	variable n, isChar
	string returnStr
	
	if (strlen(char2remove)!=1)
		Abort "The string must be one character long.  Aborting from RemoveTrailingChars"
	endif
	
	returnStr=str	
	do
		 n = strlen(returnStr)
		 if(n!=0)
		 	isChar =   ( cmpstr(returnStr[n-1], char2remove)==0  )
			if (  isChar )
				returnStr = returnStr[0, n-2]
			endif
		endif	
	while( (n!=0) &&(isChar)  )
	
	return returnStr

End


// *************************************************
// Returns a string whereby carriage return -  line feeds are removed from the end of the string

Function/s RemoveTrailingCRorLF(str)
	String str
	
	string returnStr = str
	variable n
	
	do
		n = strlen(returnStr)
		returnStr = RemoveTrailingChars("\r", returnStr)
		returnStr = RemoveTrailingChars("\n", returnStr)
	while (n != strlen(returnStr))
	
	return returnStr
	
End


// *************************************************
// Substitutes srcPat for destPat in an entire notebook.

Function ReplaceInNotebook(fileName, srcPat, destPat)
	String fileName, srcPat, destPat

	string str
	variable i = 1
	
	do 
		Notebook $fileName, selection={startOfFile, endOfFile }
	
		GetSelection notebook, $fileName,2
	
		Print "// At time "+time()+" Iteration number "+num2str(i) +" in replacing  "+srcPat+" with "+destPat
		
		str = ReplaceString(srcPat, S_Selection, destPat)  // StrSubstitute(srcPat, S_Selection, destPat)
	
		Notebook Notebook1, text = str
	
	while(i==0) //strsearch(str,srcPat,0)>=0)

	Print "// Finished at time "+time()

End


// *************************************************
// Returns a string whereby the initial chars are removed but remaining line feeds are replaced with one space.

Function/S ReplaceTheseCharsWithSpace(char2remove, str)
	String char2remove
	String str

	variable i=0, n
	string returnStr=str	
	do
		n=strlen(returnStr)
		if (n!=0)
			i = strsearch (returnStr, char2remove, i)
			if (i==0)
				returnStr = returnStr[1, n-1]
			elseif(i==(n-1))
				returnStr = returnStr[0, n-2]			
			elseif(i>0)
				returnStr = returnStr[0, i-1]+" "+str[i+1, n-1]			
			endif
		endif
	while(  (n!=0) &&(i>=0))
	
	return returnStr

End


//*************************************
//  Returns a string with all the numbers at the end removed 
// Example: StripNumericalSuffix("foo900123") returns "foo"
// Example: StripNumericalSuffix("900123") returns ""

Function/S StripNumericalSuffix(str)
	String str

	variable i, n
	string s, returnStr = ""
	
	i = strsearch(str, ".", 0)
	if (i>0)
		str = str[0,i-1]
	endif
	 n = strlen(str)
	 
	for (i = n-1; i>0;i-=1)
		s = str[i]
		if (strsearch( "0123456789", s, 0) <0)		// 	this char is not a number
			returnStr = str[0,i]
			break
		endif	
	endfor
	
	return returnStr

End	


// *************************************************
// Returns a string a substring of only numeric characters (until a non-numeric char is found)  
// this would be so easy using regular expressions, but Igor doesn't support regexps
// TeaseOutNumbers("AB990 555", 0) returns 990

Function/S TeaseOutNumbers(str, pos)
	String str
	Variable pos

	string myChar, returnStr
	variable i, n

	i = pos
	n = strlen(str)
	returnStr = ""
	
	do
		myChar = str[i]
		if(FindListItem(myChar, "-;0;1;2;3;4;5;6;7;8;9")  >=0)
			returnStr +=myChar
		endif
		i+=1 	
	while( (i<n) && ((strlen(returnStr)==0) || (FindListItem(myChar, "-;0;1;2;3;4;5;6;7;8;9")  >=0))  )

	return returnStr

End


// *************************************************
// Returns the point number where str is located in text wave; -1 if doesn't exist.

Function TextIntextWave(str, textWave)
	String str
	Wave/t textWave

	variable returnVal=-1, i=0, n = numpnts(textWave)

	for (i=0;i<n;i+=1)
		if (cmpstr(str, textWave[i])==0)
			returnVal = i
			break
		endif		
	endfor

	return returnVal

End


// *************************************************
// Returns the time value of the left most point in a time series graph.
// Be aware that this function makes some assumptions about waves and axis.

Function AxisRangeBottomLeftPt()
	
	if (strsearch (Axisinfo("", "bottom"), "CWAVE:timewave", 0)>0 )	
		WAVE twave = $"timewave"
		GetAxis/Q bottom // sets the V_min variable
		return   FindFirstGE(twave,V_min, 0)
	else
		return -1
	endif
	
End


// *************************************************
// Returns the time value of the right most point in a time series graph.
// Be aware that this function makes some assumptions about waves and axis.

Function AxisRangeBottomRightPt()

	if (strsearch (Axisinfo("", "bottom"), "CWAVE:timewave", 0)>0 )	
		WAVE twave = $"timewave"
		GetAxis/Q bottom // sets the V_min variable
		return   FindFirstGE(twave,V_max, 0)
	else
		return -1
	endif
	
End


// *************************************************
// Returns a string containing the name of the x axis on which it is plotted against.

Function/S GetXAxisStr(graphStr, traceName, instance)
	String graphStr, traceName
	Variable instance

	return  StringByKey("XAXIS", TraceInfo (graphStr, traceName, instance) )

End


// *************************************************
// Returns a string containing the name of the XWAVE it is plotted against.

Function/S GetXWaveStr(graphStr, traceName, instance)
	String graphStr, traceName
	Variable instance

	return  StringByKey("XWAVE", TraceInfo (graphStr, traceName, instance) )

End


// *************************************************
// Returns a string containing the name of the y axis on which it is plotted against.

Function/S GetYAxisStr(graphStr, traceName, instance)
	String graphStr, traceName
	Variable instance

	return  StringByKey("YAXIS", TraceInfo (graphStr, traceName, instance) )

End


// *************************************************
// For all graphs, creates a textbox with the name of the data folder containing the first  y- trace that is plotted.
// Harald's marvelous contribution
// If you want the full folder path, use GetWavesDataFolder(WaveRefIndexed("",0,1),1) instead

Function AddDataFolderTextBoxAllGraphs()

	string myGraphList,  myGraphStr, tagStr
	variable myGraphVar, myGraphNum
   
	myGraphList=WinList("*",";","WIN:1")
	myGraphNum=ItemsInList(myGraphList)
	
	for(myGraphVar=0;myGraphVar<myGraphNum;myGraphVar+=1)
       	myGraphStr = StringFromList(myGraphVar,myGraphList)
       	DoWindow/F $myGraphStr
      		tagStr=GetWavesDataFolder(WaveRefIndexed("",0,1),0)
       	TextBox/C/N=text0 tagStr		// use a default textbox name of text0  will overwrite any exiting textbox named text0
        endfor

End


// *************************************************
// Creates a more customized color legend than a default igor one.
// Stolen from Tom Ryerson - thanks Tom!

Macro Color_legend_auto(name,scalewave,color_scale,orientation,sense)					//	gets scale from WaveStats operation on chosen wave
	String name = "ted", scalewave,color_scale
	Variable orientation, sense
	Prompt name, "Pick a name for this scale:"
	Prompt scalewave, "Set the range from which wave?", popup SortList(WaveList("*",";",""), ";", 4)
	Prompt color_scale, "Choose a color scale:", popup, CTabList()
	Prompt orientation, "Whaddya want: vertical or horizontal?",popup, "vertical;horizontal"
	Prompt sense, "Pick a scale direction (up or down):",popup, "up;down"
	
	silent 1
	
	WaveStats /Q $scalewave								//	use V_min and V_max as range

	name =CleanupName(name, 0 )

	if (orientation == 1)		//	if "vertical" was picked, returning a value of 1,

		if (sense == 2)		//	if "down" was picked, returning a value of 2,
			sense = -1		//	reverse scale direction
		endif
	
		Make /O /N=(1,200) $name							//	vertical legend wave, a 1x200 matrix
		SetScale /I y, v_min,v_max,"",$name					//	sets Y scale to input range
		$name = sense*y										//	positive values of "sense" plot high values on top & vice versa
		Display; AppendImage/R $name							//	creates image plot; RH axis enabled
		ModifyImage $name ctab= {*,*,$color_scale,0}				//	use chosen color scheme
		ModifyGraph tick(right)=0,tick(bottom)=3,fSize = 12,standoff=0			//	sets up font sizes, etc.
		SetAxis right v_min, v_max								//	sets axis to desired range
		ModifyGraph noLabel(bottom)=2						//	lose the bottom axis labels
		ModifyGraph axThick(bottom)=1						//	lose the bottom axis
		ModifyGraph mirror=2									//	mirror axes, no ticks
		ModifyGraph width=20,height=250						//	tall and thin
		Label right "\\f01\\Z14"+name							//	labels axis with name of image matrix
	
	endif
	
	if (orientation == 2)		//	if "horizontal" was picked, returning a value of 2,
	
		if (sense == 2)			//	if "down" was picked, returning a value of 2,
			sense = -1			//	reverse scale direction
		endif
	
		Make /O /N=(200,1) $name							//	horizontal legend wave, a 200x1 matrix
		SetScale /I x, v_min,v_max,"",$name					//	sets bottom scale to input range
		$name = sense*x										//	positive values of "sense" plot high values on right & vice versa
		Display; AppendImage/T $name							//	creates image plot; bottom axis enabled
		ModifyImage $name ctab= {*,*,$color_scale,0}				//	use chosen color scheme
		ModifyGraph tick(top)=0,tick(left)=3,fSize = 12,standoff=0			//	sets up font sizes, etc.
		SetAxis top v_min, v_max								//	sets axis to desired range
		ModifyGraph noLabel(left)=2							//	lose the bottom axis labels
		ModifyGraph axThick(left)=1							//	lose the bottom axis
		ModifyGraph mirror=2									//	mirror axes, no ticks
		ModifyGraph width=250,height=20						//	short and wide
		Label top "\\f01\\Z14"+name							//	labels axis with name of image matrix
	
	endif
	
	//	treat the resulting picture as a normal Igor graph - doink with it as you like

End


// *************************************************
// Creates a more customized color legend than a default igor one.
// Stolen from Tom Ryerson - thanks Tom!

Macro Color_legend_manual(name,scale,v_min,v_max,orientation,sense)			
	String name = "ted",scale
	Variable v_min=0, v_max = 100, orientation, sense
	Prompt name, "Pick a name for this scale (and wave that it uses):"
	Prompt v_min, "Minimum value:"
	Prompt v_max, "Maximum value:"
	Prompt scale, "Choose a color scale:", popup, CTabList()	
	Prompt orientation, "Whaddya want: vertical or horizontal?",popup, "vertical;horizontal"
	Prompt sense, "Pick a scale direction (up or down):",popup, "up;down"

	silent 1
	
	name =CleanupName(name, 0 )	
	
	if (orientation == 1)		//	if "vertical" was picked, returning a value of 1,

		if (sense == 2)		//	if "down" was picked, returning a value of 2,
			sense = -1		//	reverse scale direction
		endif
	
		Make /O /N=(1,200) $name							//	vertical legend wave, a 1x200 matrix
		SetScale /I y, v_min,v_max,"",$name					//	sets Y scale to input range
		$name = sense*y										//	positive values of "sense" plot high values on top & vice versa
		Display; AppendImage/R $name							//	creates image plot; RH axis enabled
		ModifyImage $name ctab= {*,*,$scale,0}				//	use chosen color scheme
		ModifyGraph tick(right)=0,tick(bottom)=3,fSize = 12,standoff=0			//	sets up font sizes, etc.
		SetAxis right v_min, v_max								//	sets axis to desired range
		ModifyGraph noLabel(bottom)=2						//	lose the bottom axis labels
		ModifyGraph axThick(bottom)=1						//	lose the bottom axis
		ModifyGraph mirror=2									//	mirror axes, no ticks
		ModifyGraph width=20,height=250						//	tall and thin
		Label right "\\f01\\Z14"+name							//	labels axis with name of image matrix
	
	endif
	
	if (orientation == 2)		//	if "horizontal" was picked, returning a value of 2,
	
		if (sense == 2)			//	if "down" was picked, returning a value of 2,
			sense = -1			//	reverse scale direction
		endif
	
		Make /O /N=(200,1) $name							//	horizontal legend wave, a 200x1 matrix
		SetScale /I x, v_min,v_max,"",$name					//	sets  scale to input range
		$name = sense*x										//	positive values of "sense" plot high values on right & vice versa
		Display; AppendImage/T $name							//	creates image plot; top axis enabled
		ModifyImage $name ctab= {*,*,$scale,0}				//	use chosen color scheme
		ModifyGraph tick(top)=0,tick(left)=3,fSize = 12,standoff=0			//	sets up font sizes, etc.
		SetAxis top v_min, v_max								//	sets axis to desired range
		ModifyGraph noLabel(left)=2							//	lose the bottom axis labels
		ModifyGraph axThick(left)=1							//	lose the bottom axis
		ModifyGraph mirror=2									//	mirror axes, no ticks
		ModifyGraph width=250,height=20						//	short and wide
		Label top "\\f01\\Z14"+name							//	labels axis with name of image matrix
	
	endif
	
	//	treat the resulting picture as a normal Igor graph - doink with it as you like
	
End


// *************************************************
// Creates a more customized size legend than a default igor one.
// Stolen from Tom Ryerson - thanks Tom!

Macro Size_legend_auto(name,markernum,msize_min,msize_max,scalewave,orientation)	
	String name = "ted",scalewave
	Variable markernum,msize_min=1, msize_max=8,orientation
	Prompt name, "Pick a name for this scale:"
	Prompt markernum, "Choose a marker:", popup, "+;x;circle;square;triangle;filled circle;filled square;filled triangle"
	Prompt msize_min, "Smallest marker size:"
	Prompt msize_max, "Largest marker size:"
	Prompt scalewave, "Set the range from which wave?", popup SortList(WaveList("*",";",""), ";", 4)
	Prompt orientation, "Whaddya want: vertical or horizontal?",popup, "vertical;horizontal"
	
	silent 1;PauseUpdate
	
	WaveStats /Q $scalewave	// use v_min and v_max from this wave to set range scale
	
	variable markertype = 0
	string Xname, Yname
	if (markernum == 1)		//	"+" selected, returning a value of 1,
		markertype = 0				//	chooses "+" marker (see Ygor manual, p. II-301, for table of marker numbers).
	endif
	if (markernum == 2)		//	"x" selected
		markertype = 1
	endif
	if (markernum == 3)		//	"circle" selected
		markertype = 8
	endif
	if (markernum == 4)		//	"square" selected
		markertype = 5
	endif
	if (markernum == 5)		//	"triangle" selected
		markertype = 6
	endif
	if (markernum == 6)		//	"filled circle" selected
		markertype = 19
	endif
	if (markernum == 7)		//	"filled square" selected
		markertype = 16
	endif
	if (markernum == 8)		//	"filled triangle" selected
		markertype = 17
	endif
	
	Yname = CleanupName(name+"_Yscale",0)
	Xname = CleanupName(name+"_Xscale",0)

	if (orientation == 1)		//	if "vertical" was picked, returning a value of 1,
		
		Make/o/d/n=10/D $Yname,$Xname												//	10-point waves
		$Yname = (v_min + ((v_max-v_min)/20)) + (p * (v_max-v_min)/10)		//	gets equally-spaced points throughout Zscale_wave range
		$Xname = 0																		//	can be any value, really -
			
		Display/R $Yname vs $Xname					//	plots on RH axis
		ModifyGraph mode=3							//	plot as markers
		ModifyGraph marker=markertype				//	marker selection from input above
		ModifyGraph rgb=(0,0,0)						//	black color
		
			//	zmrkSize={zWave,zMin,zMax,mrkmin,mrkmax }
		ModifyGraph zmrkSize($Yname)={$Yname,v_min,v_max,msize_min,msize_max}	//	"*" autoscales based on $Xname range
		
		ModifyGraph width=20,height=250				//	short and wide
		ModifyGraph mirror=2							//	frames everything nicely
		ModifyGraph fSize=12,standoff=0				//	12-pt axis labels, no standoff
		SetAxis right v_min, v_max						//	sets to desired axis range
		ModifyGraph tick(bottom)=3					//	no bottom axis ticks
		ModifyGraph noLabel(bottom)=2				//	no bottom axis tick labels
		Label right "\\f01\\Z14"+name					//	labels axis with name of legend
	
	endif	
	
	if (orientation == 2)		//	if "horizontal" was picked, returning a value of 2,
		
		Make/o/d/n=10/D $Yname,$Xname												//	10-point waves
		$Yname = 0																		//	can be any value, really -
		$Xname = (v_min + ((v_max-v_min)/20)) + (p * (v_max-v_min)/10)		//	gets equally-spaced points throughout range
		
		Display /T $Yname vs $Xname					//	displays on top axis
		ModifyGraph mode=3							//	plot as markers
		ModifyGraph marker=markertype				//	marker selection from input above
		ModifyGraph rgb=(0,0,0)						//	black color
		
			//	zmrkSize={zWave,zMin,zMax,mrkmin,mrkmax }
		ModifyGraph zmrkSize($Yname)={$Xname,v_min,v_max,msize_min,msize_max}	//	"*" autoscales based on $Xname range
		
		ModifyGraph width=250,height=20				//	short and wide
		ModifyGraph mirror=2							//	frames everything nicely
		ModifyGraph fSize=12,standoff=0				//	12-pt axis labels, no standoff
		SetAxis top v_min, v_max						//	sets to desired axis range
		ModifyGraph tick(left)=3						//	no LH axis ticks
		ModifyGraph noLabel(left)=2					//	no LH axis tick labels
		Label top "\\f01\\Z14"+name					//	labels axis with name of legend
		
	endif
	
End


// *************************************************
// Creates a more customized size legend than a default igor one.
// Stolen from Tom Ryerson - thanks Tom!

Macro Size_legend_manual(name,markernum,msize_min,msize_max,scale_min,scale_max,orientation)	
	String name = "ted"
	Variable markernum,msize_min=0.25, msize_max=8,scale_min=0,scale_max=100,orientation
	Prompt name, "Enter name for scale (29 chrs max):"
	Prompt markernum, "Choose a marker:", popup, "+;x;circle;square;triangle;filled circle;filled square;filled triangle"
	Prompt msize_min, "Smallest marker size:"
	Prompt msize_max, "Largest marker size:"
	Prompt scale_min, "Minimum axis scale:"
	Prompt scale_max, "Maximum axis scale:"
	Prompt orientation, "Whaddya want: vertical or horizontal?",popup, "vertical;horizontal"
	
	silent 1;PauseUpdate
	
	variable markertype, scale_range
	scale_range = scale_max - scale_min
	
	string Xname, Yname 
	
	if (strlen(name)>29)
		Abort "String entered for name of scale is too long. 29 characters is the maximum."
	endif

	if (markernum == 1)		//	"+" selected, returning a value of 1,
		markertype = 0				//	chooses "+" marker (see manual, p. II-301, for table of marker numbers).
	endif
	
	if (markernum == 2)		//	"x" selected
		markertype = 1
	endif
		
	if (markernum == 3)		//	"circle" selected
		markertype = 8
	endif
		
	if (markernum == 4)		//	"square" selected
		markertype = 5
	endif
	
	if (markernum == 5)		//	"triangle" selected
		markertype = 6
	endif
		
	if (markernum == 6)		//	"filled circle" selected
		markertype = 19
	endif
	
	if (markernum == 7)		//	"filled square" selected
		markertype = 16
	endif
	
	if (markernum == 8)		//	"filled triangle" selected
		markertype = 17	
	endif
	
	Yname = CleanupName(name+"_Y",0)
	Xname =  CleanupName(name+"_X",0)
	
	if (orientation == 1)		//	if "vertical" was picked, returning a value of 1,
	
		Make/o/d/n=10/D $Yname,$Xname										//	10-point waves
		$Yname = (scale_min + scale_range/20) + (p * (scale_range/10))	//	gets equally-spaced points throughout range
		$Xname = 0																//	can be any value, really -
			
		Display/R $Yname vs $Xname					//	plots on RH axis
		ModifyGraph mode=3							//	plot as markers
		ModifyGraph marker=markertype				//	marker selection from input above
		ModifyGraph rgb=(0,0,0)						//	black color
		
			//	zmrkSize={zWave,zMin,zMax,mrkmin,mrkmax }
		ModifyGraph zmrkSize($Yname)={$Yname,scale_min,scale_max,msize_min,msize_max}	//	"*" autoscales based on $Xname range
		
		ModifyGraph width=20,height=250				//	short and wide
		ModifyGraph mirror=2							//	frames everything nicely
		ModifyGraph fSize=12,standoff=0				//	12-pt axis labels, no standoff
		SetAxis right scale_min, scale_max				//	sets to desired axis range
		ModifyGraph tick(bottom)=3					//	no bottom axis ticks
		ModifyGraph noLabel(bottom)=2				//	no bottom axis tick labels
		Label right "\\f01\\Z14"+name					//	labels axis with name of legend

	endif
			
	if (orientation == 2)		//	if "horizontal" was picked, returning a value of 2,
				
		Make/o/d/n=10/D $Yname,$Xname										//	10-point waves
		$Yname = 0																//	can be any value, really -
		$Xname = (scale_min + scale_range/20) + (p * (scale_range/10))	//	gets equally-spaced points throughout range
		
		Display /T $Yname vs $Xname					//	displays on top axis
		ModifyGraph mode=3							//	plot as markers
		ModifyGraph marker=markertype				//	marker selection from input above
		ModifyGraph rgb=(0,0,0)						//	grey color
		
			//	zmrkSize={zWave,zMin,zMax,mrkmin,mrkmax }
		ModifyGraph zmrkSize($Yname)={$Xname,scale_min,scale_max,msize_min,msize_max}	//	"*" autoscales based on $Xname range
		
		ModifyGraph width=250,height=20				//	short and wide
		ModifyGraph mirror=2							//	frames everything nicely
		ModifyGraph fSize=12,standoff=0				//	12-pt axis labels, no standoff
		SetAxis top scale_min, scale_max				//	sets to desired axis range
		ModifyGraph tick(left)=3						//	no LH axis ticks
		ModifyGraph noLabel(left)=2					//	no LH axis tick labels
		Label top "\\f01\\Z14"+name					//	labels axis with name of legend
		
	endif
	
End


// *************************************************
// Returns a "0" if the cursors on the top graph are on the same wave  and returns a "1" if not.

Function Csr_on_wave()
	
	Variable returnVal = 0

	string waveAStr=CsrWave(A), waveBStr=CsrWave(B)
	
	if(strlen(waveAStr)==0 )
		returnVal = 1
	endif
	if (strlen(waveBStr)==0)
		returnVal = 1
	endif	
	if (cmpstr(waveAStr, waveBStr) )
		returnVal = 1
	endif	
	
	return(returnVal)

End


// *************************************************
//  Calculates an average between two cursor positions on the top graph and inserts average into a wave.
// It averages the data values of the wave in the top graph which has the string
// "counts" in it, i.e., no_counts.  It changes the wave in which the
// cursors are on, by inserting the average between the cursors.
// It also inserts nans elsewhere between the cursors and a delta-many beyond the cursors.

Macro Insert_new_avg()

	silent 1
	variable r,s, delta = 20
	
	string waveA=CsrWave(A), waveB=CsrWave(B)
	string thewavelist= WaveList("*count*",";","WIN:")
	string the_avg_wavestring=StringFromList(0, thewavelist, ";")
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "The cursors need to be on the same wave - Aborting from Insert_new_avg"
	endif
	
	if (  cmpstr((waveA[strlen(waveA)-4,strlen(waveA) ]), "_avg")    )
		Abort "Please put your cursor on a wave that ends in _avg - Aborting from Insert_new_avg"
	endif
	
	r=pcsr(A)
	s=pcsr(B)
	
	WaveStats/Q/R=[r,s] $the_avg_wavestring
	
	$waveA[r-delta,s+delta]=nan
	$waveA[(r+s)/2]=V_avg
	
	Print "// You inserted the average ",V_avg,"at point ",(r+s)/2, " in wave "+waveA

End


// *************************************************
// Insert the average values between two cursors on the top graph.
//  Overwrites any non-nan values between cursors!

Function Make_avg_between_cursors()

	variable x0, x1
	string waveAStr=CsrWave(A)
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "The cursors need to be on the same wave - Aborting from Make_avg_between_cursors"
	endif
	
	WAVE waveA = CsrWaveRef(A)

	x0=pcsr(A)
	x1=pcsr(B)
		
	if (x0 == x1 )
		Abort "The cursors need to be on different values - Aborting from Make_avg_between_cursors"
	endif
	if (x0>x1)
		abort "Cursor B needs to be ahead of cursor A."
	endif

	Wavestats/q/r=[x0, x1] waveA
	
	waveA[x0, x1]= V_avg

	Print "// You inserted averaged values of "+num2str(V_avg)+" in points ",x0," through ",x1," in wave "+waveAStr

End


// *************************************************
// Inserts linearly interpolated values between two cursors on the top graph.
//  Overwrites any non-nan values between cursors!

Function Make_interp_between_cursors()

	variable x0, x1, y0, y1, m
	string waveAStr=CsrWave(A)
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "the cursors need to be on the same wave - Aborting from Make_interp_between_cursors"
	endif
	
	WAVE waveA = CsrWaveRef(A)

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
	print m

	Print "// You inserted interpolated values in points ",x0," through ",x1," in wave "+waveAStr

End


// *************************************************
// Inserts nans between two cursors on the top graph.

Function Make_nans_between_cursors()	

	string waveAStr = CsrWave(A)
	print waveAStr
	variable r,s
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "The cursors need to be on the same wave - Aborting from Make_nans_between_cursors"
	endif
	
	WAVE waveA = CsrWaveRef(A)
	
	r=pcsr(A)
	s=pcsr(B)
	if (r>s)
		abort "Cursor B needs to be ahead of cursor A."
	endif
	
	waveA[r,s]=nan
	Print "// You naned points ",r," through ",s," in wave "+waveAStr

End


// *************************************************
// Inserts ones between two cursors on the top graph.

Function Make_ones_between_cursors()	

	string waveAStr = CsrWave(A)
	variable r,s
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "The cursors need to be on the same wave - Aborting from Make_ones_between_cursors"
	endif
	
	WAVE waveA = CsrWaveRef(A)
	
	r=pcsr(A)
	s=pcsr(B)
	if (r>s)
		abort "Cursor B needs to be ahead of cursor A."
	endif
	
	waveA[r,s]=1
	Print "// You inserted 1s in points ",r," through ",s," in wave "+waveAStr

End


// *************************************************
// Inserts zeros between two cursors on the top graph.

Function Make_zeros_between_cursors()

	string waveAStr = CsrWave(A)
	variable r,s
	
	if( csr_on_wave()  )		//  csr_on_wave()  returns 1 if you need to abort
		Abort "the cursors need to be on the same wave - Aborting from Make_zeros_between_cursors"
	endif
	
	WAVE waveA = CsrWaveRef(A)
	
	r=pcsr(A)
	s=pcsr(B)
	if (r>s)
		abort "Cursor B needs to be ahead of cursor A."
	endif
	
	waveA[r,s]=0
	Print "// You inserted zeros in points ",r," through ",s," in wave "+waveAStr

End


// *************************************************
// Inserts nan values into wave yWaveStr and its matching XWave based on marquee settings.
// Users should be aware that this function makes assumptions about axis names.
// Assumes that yWave is plotted against the left and X is plotted against the bottom axes
// If one of the two waves has nan values the other will not be naned, even though it's value may be within the marquee range.

Function NanInsideMarqueeXYData(yWaveStr)
	String yWaveStr

	Wave yWave = TraceNameToWaveRef("", yWaveStr )	// assumes the top graph. 
	if (WaveExists(yWave)==0)
		Abort "The wave "+yWaveStr+"does not exist or is not on the graph. - Aborting from NanInsideMarqueeXYData"
	endif
	
	variable numNaned=0, i, n=numpnts(yWave)
	
	WAVE/Z xWave = XWaveRefFromTrace("", yWaveStr)
	if (WaveExists(xWave)==0)
		Abort "The wave is not plotted against an x wave. - Aborting from NanInsideMarqueeXYData"
	endif
	
	GetMarquee left, bottom			// find vertical marquee location
	
	for (i=0;i<n;i+=1)
		if  (   (yWave[i]<V_top)&& (yWave[i]>V_bottom)&& (xWave[i]<V_right)&&(xWave[i]>V_left) )
			numNaned+=1
			yWave[i]=nan
			xWave[i]=nan
		endif
	endfor			
	
	Print "// You have naned "+num2str(numNaned)+ " points between "+num2str(V_bottom)+" and  "+num2str(V_top)+" in wave "+yWaveStr

End


// *************************************************
// Inserts nan values into wave yWaveStr based on marquee settings.
// Users should be aware that this function makes assumptions about axis names.

Function NanInsideMarqueeYData(yWaveStr)
	String yWaveStr

	Wave yWave = TraceNameToWaveRef("", yWaveStr )	// assumes the top graph. 
	
	if (WaveExists(yWave)==0)
		Abort "The wave "+yWaveStr+"does not exist or is not on the graph - Aborting from NanInsideMarqueeYData"
	endif
	
	variable numNaned=0, i, n=numpnts(yWave)
	variable minPoint, maxPoint
	variable minY, maxY
	variable isParametricPlot
	
	// find the wave against which yWave is plotted, if any
	WAVE/Z xWave = XWaveRefFromTrace("", yWaveStr)
	isParametricPlot = WaveExists(xWave)	// 1 if exists, 0 if not
	
	GetMarquee left				// find vertical marquee location is in terms of the left axis
	minY = V_bottom
	maxY = V_top
	
	GetMarquee bottom					// find horizontal marquee location is in terms of the bottom axis
	
	if (isParametricPlot)
		for (i=0;i<n;i+=1)
			if  (   (yWave[i] <V_top)&& (yWave[i] > V_bottom)&& (xWave[i] <V_right)&&(xWave[i] >V_left) )
				numNaned+=1
				yWave[i]=nan
			endif			
		endfor
		Print "// You have naned "+num2str(numNaned)+ " points between "+num2str(V_bottom)+" and  "+num2str(V_top)+" in wave "+yWaveStr
	else		// faster this way
		minPoint = ceil(x2pnt(yWave, V_left))
		maxPoint = floor(x2pnt(yWave, V_right))
		yWave[minPoint, maxPoint]=nan
		Print "// You have naned the points "+num2str(minPoint)+ " through "+num2str(maxPoint)+" in wave "+yWaveStr
	endif

End


// *************************************************
// Makes an arbitrary number of different colors in a color table.
// Be aware that all the colors may not be distinct at a distance.

Function MakeRainbowColorMany(n)
	Variable n

	variable maxVal, index, num
	
	num = n + (27 - mod(n, 27))
	maxVal = 65535
	
	Make/o/d/n=(num,3) RainbowColorTableMany

	for (index = 1;index <= num/27; index+=1)
		
		RainbowColorTableMany[0+(index-1)*27]={ {maxVal/index}, {0}, {0} }		//red
		RainbowColorTableMany[1+(index-1)*27]={ {0}, {maxVal/index}, {0} }		//green
		RainbowColorTableMany[2+(index-1)*27]={ {0}, {0}, {maxVal/index} }		//blue
	
		RainbowColorTableMany[3+(index-1)*27]={ {maxVal/index}, {maxVal/index}, {0} }	//yellow
		RainbowColorTableMany[4+(index-1)*27]={ {0}, {maxVal/index}, {maxVal/index} }	//cyan
		RainbowColorTableMany[5+(index-1)*27]={ {maxVal/index}, {0}, {maxVal/index} }	//magenta
	
		RainbowColorTableMany[6+(index-1)*27]={ {maxVal/(index*2)}, {0}, {0} }	//brick
		RainbowColorTableMany[7+(index-1)*27]={ {0}, {maxVal/(index*2)}, {0} }	//darkgreen
		RainbowColorTableMany[8+(index-1)*27]={ {0}, {0}, {maxVal/(index*2)} }	//darkblue
	
		RainbowColorTableMany[9+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/(index*2)}, {0} }	//mud
		RainbowColorTableMany[10+(index-1)*27]={ {0}, {maxVal/(index*2)}, {maxVal/(index*2)} }	//teal
		RainbowColorTableMany[11+(index-1)*27]={ {maxVal/(index*2)}, {0}, {maxVal/(index*2)} }	//redpurple
	
		RainbowColorTableMany[12+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/index}, {0} }	//lightgreen
		RainbowColorTableMany[13+(index-1)*27]={ {0}, {maxVal/(index*2)}, {maxVal/index} }	//medblue
		RainbowColorTableMany[14+(index-1)*27]={ {maxVal/index}, {0}, {maxVal/(index*2)} }	//pinkyred
	
		RainbowColorTableMany[15+(index-1)*27]={ {maxVal/index}, {maxVal/(index*2)}, {0} }	//orange
		RainbowColorTableMany[16+(index-1)*27]={ {0}, {maxVal/index}, {maxVal/(index*2)} }	//othergreen
		RainbowColorTableMany[17+(index-1)*27]={ {maxVal/(index*2)}, {0}, {maxVal/index} }	//darkpurple
	
		RainbowColorTableMany[18+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/index}, {maxVal/(index*2)} }	//palegreen
		RainbowColorTableMany[19+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/(index*2)}, {maxVal/index} }	//slate
		RainbowColorTableMany[20+(index-1)*27]={ {maxVal/index}, {maxVal/(index*2)}, {maxVal/(index*2)} }	//coral
	
		RainbowColorTableMany[21+(index-1)*27]={ {maxVal/index}, {maxVal/(index*2)}, {maxVal/index} }	//pink
		RainbowColorTableMany[22+(index-1)*27]={ {maxVal/index}, {maxVal/index}, {maxVal/(index*2)} }	//lightyellow
		RainbowColorTableMany[23+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/index}, {maxVal/index} }	//paleblue
	
		RainbowColorTableMany[24+(index-1)*27]={ {0}, {0}, {0} }					//black
		RainbowColorTableMany[25+(index-1)*27]={ {maxVal/(index*2)}, {maxVal/(index*2)}, {maxVal/(index*2)} }		//gray
		RainbowColorTableMany[26+(index-1)*27]={ {maxVal/index}, {maxVal/index}, {maxVal/index} }		//white

	endfor

End


// *************************************************
// Makes a color table of 21 different colors

Function MakeRainbowColorTable()

	Make/o/d/n=(21,3) RainbowColorTable

	RainbowColorTable[0]={ {65535}, {0}, {0} }		//red
	RainbowColorTable[1]={ {0}, {65535}, {0} }		//green
	RainbowColorTable[2]={ {0}, {0}, {65535} }		//blue

	RainbowColorTable[3]={ {65535}, {65535}, {0} }	//yellow
	RainbowColorTable[4]={ {0}, {65535}, {65535} }	//cyan
	RainbowColorTable[5]={ {65535}, {0}, {65535} }	//magenta

	RainbowColorTable[6]={ {32767}, {0}, {0} }	//brick
	RainbowColorTable[7]={ {0}, {32767}, {0} }	//darkgreen
	RainbowColorTable[8]={ {0}, {0}, {32767} }	//darkblue

	RainbowColorTable[9]={ {32767}, {32767}, {0} }	//mud
	RainbowColorTable[10]={ {0}, {32767}, {32767} }	//teal
	RainbowColorTable[11]={ {32767}, {0}, {32767} }	//redpurple

	RainbowColorTable[12]={ {32767}, {65535}, {0} }	//lightgreen
	RainbowColorTable[13]={ {0}, {32767}, {65535} }	//medblue
	RainbowColorTable[14]={ {65535}, {0}, {32767} }	//pinkyred

	RainbowColorTable[15]={ {65535}, {32767}, {0} }	//orange
	RainbowColorTable[16]={ {0}, {65535}, {32767} }	//othergreen
	RainbowColorTable[17]={ {32767}, {0}, {65535} }	//darkpurple

	RainbowColorTable[18]={ {32767}, {65535}, {32767} }	//palegreen
	RainbowColorTable[19]={ {32767}, {32767}, {65535} }	//slate
	RainbowColorTable[20]={ {65535}, {32767}, {32767} }	//coral

	RainbowColorTable[21]={ {65535}, {32767}, {65535} }	//pink
	RainbowColorTable[22]={ {65535}, {65535}, {32767} }	//lightyellow
	RainbowColorTable[23]={ {32767}, {65535}, {65535} }	//paleblue

	RainbowColorTable[24]={ {0}, {0}, {0} }					//black
	RainbowColorTable[25]={ {32767}, {32767}, {32767} }		//gray
	RainbowColorTable[26]={ {65535}, {65535}, {65535} }		//white

End


// *************************************************
// Similar to Rainbowize traces, except that it will rainbowize an arbitrarily large number of traces.

Function RainbowizeManyTraces()

	DoUpdate
	
	string traceStr, traceList
	variable i,n 

	traceList = TraceNameList("", ";", 1)
	n = ItemsInList(traceList)
	MakeRainbowColorMany(n)
	wave RainbowColorTableMany = $"RainbowColorTableMany"
	
	for (i=0; i<n; i+=1)
		traceStr = StringFromList(i,traceList,  ";")
		ModifyGraph rgb($traceStr)=(RainbowColorTableMany[i][0],RainbowColorTableMany[i][1],RainbowColorTableMany[i][2]);DelayUpdate
	endfor
	
	KillwavesinList("RainbowColorTableMany")

End


// *************************************************
// Will change the colors of the first 27 traces on the topmost graph

Function RainbowizeTraces()

	DoUpdate
	
	string traceStr, traceList = TraceNameList("", ";", 1)
	variable i,n = ItemsInList(traceList)

	MakeRainbowColorTable()
	WAVE RainbowColorTable = $"RainbowColorTable"
	
	for (i=0; i<min(n, 26); i+=1)
		traceStr = StringFromList(i,traceList,  ";")
		ModifyGraph rgb($traceStr)=(RainbowColorTable[i][0],RainbowColorTable[i][1],RainbowColorTable[i][2]);DelayUpdate
	endfor
	
	KillwavesinList("RainbowColorTable")

End


// *************************************************
// A specialized function used for appending time tags to lat -lon plots.
// It makes assumptions about the names of waves and that the timewave is in 1 Hz.

Function/S AppendTimeTags2LatLong(timeWaveStr, gpsLatWaveStr, minuteInterval)
	String timeWaveStr, gpsLatWaveStr
	Variable minuteInterval

	if (  (Exists(timeWaveStr)!=1) || (Exists(gpsLatWaveStr)!=1) || (minuteInterval <=0)  )
		minuteInterval = 60
		prompt timeWaveStr,  "The timewave",popup,WaveList("*time*",";","")
		prompt gpsLatWaveStr,  "The latwave",popup,WaveList("*at*",";","WIN:")
		prompt minuteInterval "The number of minutes as an interval for the tags"
		doPrompt "Enter the values.", timeWaveStr, gpsLatWaveStr, minuteInterval
	endif

	WAVE twave = $timeWaveStr
	WAVE gpsLatWave = $gpsLatWaveStr

	Duplicate/o twave timeKey
	
	timeKey = SelectNumber(mod(twave[p], 60*minuteInterval)==0, nan, twave[p] ) 

	variable k=0, i=0, n=numpnts(twave)
	string nameStr

	k = FindFirst(timeKey, k+1)

	do
//		print k
		nameStr = "time"+num2str(i)
		Tag/I=1/C/N=$nameStr/X=4/Y=4/A=MC/L=1/F=0/M $gpsLatWaveStr, k, "\\F'Times'"+Secs2Time(twave[k], 2)
		k = FindFirst(timeKey, k+1)	
		i+=1	
	while(k<n)
	
	Killwaves timeKey
	
End


// *************************************************
// A specialized function used for appending time tags and a textbox with a date to lat -lon plots.
// It makes assumptions about the names of waves and that the timewave is in 1 Hz.

Function AppendTimeTags2LatLongWithDate(timeWaveStr, gpsLatWaveStr, minuteInterval)
	String timeWaveStr, gpsLatWaveStr
	Variable minuteInterval
	
	AppendTimeTags2LatLong(timeWaveStr, gpsLatWaveStr, minuteInterval)
	
	WAVE twave = $timeWaveStr
	
	variable startTime = twave[0]
	string dateStr = time2yymmddhhmmss(startTime)
	
	dateStr = dateStr[0,5]
	
	TextBox/A=RT/C/N=fltdateBox/S=0 "\Z14flt"+dateStr

End

// *************************************************
// Wind conversion macros.
// *************************************************
// One should convert traditional wind direction and speed
// data to wind north/east component vectors before doing any averaging.
// The averaged data can be converted back into the
//  traditional wind speed and direction vectors after averaging.

// Averages wind direction wave by first converting wind direction and speed into north and 
// east components, then averaging (currently only avgxsec is available)
// and then converting averaged components back into true averaged wind direction.
// This way, when wind direction is from near north, no arbitrary average 
// southerly winds are produced, as they would have been if wind direction would have
// been averaged directly.
// Calls Wind_spd_dir_2_n_e() and Wind_n_e_2_dir()

Macro avg_WindDir_Xsec(WindDirStr, WindSpdStr, TimeWaveStr, avgXsec)
	String WindDirStr = "WindDir"
	String WindSpdStr = "WindSpd"
	String TimeWaveStr = "AOCTimewave"
	Variable avgXsec = 60
	Prompt WindDirStr,"Wind Direction wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt WindSpdStr,"Wind speed wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt TimeWaveStr,"Time wave",popup SortList(WaveList("*",";",""), ";", 16)

	// convert wind direction and speed to north and east components
	Wind_spd_dir_2_n_e("w_e","w_n",windDirStr,windSpdStr)
	String WindAvgStr = CleanupName(WindDirStr+"_avg",0)
	// average north and east components
	AveragOnlyXSecs(timeWaveStr, "w_e", "w_e_avg", avgXsec)
	AveragOnlyXSecs(timeWaveStr, "w_n", "w_n_avg", avgXsec)
	
	//Convert averaged components to wind direction
	Wind_n_e_2_dir("w_e_avg","w_n_avg",WindAvgStr)
	
	// Make mid time wave
	String StartTimeStr = "start_"+num2str(avgXsec)+"sec"
	String StopTimeStr = "stop_"+num2str(avgXsec)+"sec"
	String MidTimeStr = "mid_"+num2str(avgXsec)+"sec"
	Duplicate/O $StartTimeStr, $MidTimeStr
	$MidTimeStr = (($StartTimeStr) + ($StopTimeStr))/2

End Macro

// *************************************************
// Converts wind north and east component vectors (waves)
// into wind speed and wind direction waves.
// Prompts the user for the names for the newly created
// direction and speed vectors... has default names of "w_dir" and "w_spd"
// Called by AddWindbarbs.

Macro Wind_n_e_2_spd_dir(windEastStr,windNorthStr,windDirStr,windSpdStr)
	String windEastStr="w_e",windNorthStr="w_n",windDirStr="w_dir",windSpdStr="w_spd"
	Prompt windEastStr,"East component wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windNorthStr,"North component wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windDirStr,"Name of direction wave"
	Prompt windSpdStr,"Name of velocity wave"

	windDirStr = CleanupName(windDirStr, 0)
	windSpdStr = CleanupName(windSpdStr, 0)
	
	Duplicate/o $windEastStr, $windDirStr, $windSpdStr	
	
	$windDirStr=mod(atan2($windEastStr,$windNorthStr)/pi*180+360,360)
	$windSpdStr=sqrt( ($windNorthStr*$windNorthStr)+($windEastStr*$windEastStr))

End

//	Called by avg_WindDir_Xsec()
Macro Wind_n_e_2_dir(windEastStr,windNorthStr,windDirStr)
	String windEastStr="w_e",windNorthStr="w_n",windDirStr="w_dir"
	Prompt windEastStr,"East component wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windNorthStr,"North component wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windDirStr,"Name of direction wave"

	windDirStr = CleanupName(windDirStr, 0)
	
	Duplicate/o $windEastStr, $windDirStr	
	
	$windDirStr=mod(atan2($windEastStr,$windNorthStr)/pi*180+360,360)

End

// *************************************************
// Converts wind speed and wind direction waves
// into wind north and east component vectors.
// Prompts the user for the names for the newly created
// north and east vectors... has default names of w_e" and "w_n"

Macro Wind_spd_dir_2_n_e(windEastStr,windNorthStr,windDirStr,windSpdStr)
	String windEastStr="w_e",windNorthStr="w_n",windDirStr="w_dir",windSpdStr="w_spd"
	Prompt windDirStr,"Wind Direction wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windSpdStr,"Wind velocity wave",popup SortList(WaveList("*",";",""), ";", 16)
	Prompt windEastStr,"Name to give east component wave"
	Prompt windNorthStr,"Name to give north component wave"

	windNorthStr = CleanupName(windNorthStr, 0)
	windEastStr = CleanupName(windEastStr, 0)
	
	Duplicate/o $windDirStr, $windEastStr, $windNorthStr
	
	$windNorthStr=cos($windDirStr/180*pi)*$windSpdStr
	$windEastStr=sin($windDirStr/180*pi)*$windSpdStr

End


// *************************************************
//  Contributed by Tara Fortin, edited by Donna Sueper
//
// This macro creates wind barbs for 60 second averages and adds them to map graph.
// If something other than 60 second averages is needed, need to change numrows and k and l.
// If no graph called "Map" exists, need to change appendtograph line.
// Scale is such that flags match convention (see e.g. Wallace and Hobbs, p.113). See notes below.
//
// Some things to note. Conventional meteorological flags are in knots and change every 5 knots.
// Igor goes from 0 to 40 and changes every 1. So Igor's flag at 1 is the same as the conventional flag
// at 5 knots. So I convert averaged wind speed to knots and then divide by 5.
// Also, met data for direction is 0 to 360 degrees with 0 corresponding to N and moving clockwise. 
// Igor wants things in radians and has zero where met data has 90 degrees. Also, Igor moves
// counterclockwise. So to get all markers facing the correct direction, I took the average direction,
// subtracted 90 degrees to adjust for the offset in the zero position, and flipped by multiplying by -1.
//
//  Does not rely on the data being in 1 sec time resolution
//  Creates or overwrites the following waves: windBarbLat, windBarbLon, windBarbMatrix, windDirAvg, windSpdAvg.

Macro AddWindBarbs(windDirStr,windSpdStr,latStr,lonStr, timeStr, XSecs2avg)
	String windDirStr,windSpdStr,latStr,lonStr, timeStr
	Variable  XSecs2avg=60		// a default
	prompt windDirStr,"The wind direction wave?", popup SortList(WaveList("*",";",""), ";", 16)// typically windDir is used; can be replaced with WaveList("wind*",";","")
	prompt windSpdStr,"The wind speed wave?", popup SortList(WaveList("*",";",""), ";", 16)// typically windSpd is used; can be replaced with WaveList("wind*",";","")
	prompt latStr,"The latitude wave?", popup SortList(WaveList("*",";",""), ";", 16)		// typically gpsLat is used; can be replaced with WaveList("G*",";","")
	prompt lonStr,"The longitude wave?", popup SortList(WaveList("*",";",""), ";", 16)		// typically gpsLon is used; can be replaced with WaveList("G*",";","")
	prompt timeStr,"The time wave?", popup SortList(WaveList("*",";",""), ";", 16)// typically windDir is used; can be replaced with WaveList("*",";","")
	prompt XSecs2avg, "Secs to avg and display barbs?"
	
	silent 1;PauseUpdate
	
	variable delta, n
	
	Print "AddWindBarbs(\""+windDirStr+"\",\""+windSpdStr+"\",\""+latStr+"\",\""+lonStr+"\",\""+timeStr+"\","+num2str(XSecs2avg)+")"

	if ( (Exists(windDirStr)!=1) || (Exists(windSpdStr)!=1) || (Exists(latStr)!=1) || (Exists(lonStr)!=1) || (Exists(timeStr)!=1)  )
		Abort "One of these waves do not exist: "+ windDirStr+windSpdStr+latStr+lonStr+ timeStr+ "- Aborting from AddWindBarbs"
	endif
	
	Wind_spd_dir_2_n_e("windEast","windNorth",windDirStr,windSpdStr)		// a call to another macro; assumes these default waves names are ok
	
	AveragOnlyXSecs(timeStr, "windEast", "windEastAvg", XSecs2avg)
	AveragOnlyXSecs(timeStr, "windNorth", "windNorthAvg", XSecs2avg)
	
	Wind_n_e_2_spd_dir("windEastAvg","windNorthAvg","windDirAvg","windSpdAvg") 	// a call to another macro; assumes these default waves names are ok

	Killwaves/Z windEast,windNorth			// one can comment out if preferred
	Killwaves/Z windEastAvg, windNorthAvg 	// one can comment out if preferred

	n = numpnts(windDirAvg)
	if (n*XSecs2avg > numpnts($latstr))	//Can't exceed number of points in $latstr and $lonstr
		n-=1
	endif
	
	Make/D/O/N=(n,3) windBarbMatrix
	Make/D/O/N=(n) windBarbLat,windBarbLon

	windBarbMatrix[][0] =  20
	windBarbMatrix[][1] = ((-1*(windDirAvg[p]-90))*pi)/180  // convert to radians
	windBarbMatrix[][2] = ((windSpdAvg[p])*1.95)/5	 //convert to knots and then igor windbarb scale
	SetDimLabel 1,2,windBarb,windBarbMatrix
	
	Make/o/d/n=(n) timeDiffTemporaryWave
	timeDiffTemporaryWave[1, ] = ($timeStr[p] - $timeStr[p-1])	//start on point 1, not 0 since p-1 was out of range.
	timeDiffTemporaryWave[0] = nan
	WaveStats/Q timeDiffTemporaryWave
	Killwaves timeDiffTemporaryWave
	
	if (XSecs2avg==1)								// Then just get every lat/lon pair. Fixes problem when XSecs2avg=1
		windBarbLat = $latStr
		windBarbLon = $lonStr
	else
		if (V_min == V_max)						// time is in even increments, pick out a middle value
			delta = round( XSecs2avg/V_min )	
			windBarbLat = $latStr[round(delta/2) + p*delta]
			windBarbLon = $lonStr[round(delta/2) + p*delta]
		else											// time is not in even increments, average lat and lon
			AveragOnlyXSecs(timeStr, latStr, "windBarbLat", XSecs2avg)
			AveragOnlyXSecs(timeStr, lonStr, "windBarbLon", XSecs2avg)
		endif
	endif
	RemoveFromGraph/Z windBarbLat                        // Remove it if it is already on the graph 
	AppendToGraph windBarbLat vs windBarbLon			// append to top graph
	ModifyGraph mode(windBarbLat)=3,arrowMarker(windBarbLat)={windBarbMatrix,1,5,0.5,0}
	ModifyGraph rgb(windBarbLat)=(43690,43690,43690)		// make barbs default color gray

End

//	Add 2 functions to marquee menu
Menu "GraphMarquee"
	"--"
	"SaveAxisSettings"
	"RestoreAxisSettings"
End

//	Run this function to save the axis min & max for all axes on graph (no matter how many there are).
//	Then Expand, Shrink, pan, etc. as user wishes
//	To get back to original graph before Expand, Shrink, etc. run "restoreAxisSettings()"
//	This gives the functionality of multiple undos.
//Makes string called gAxisSettingsStr - formatted like: axisname0;min0;max0;axisname1;min1;max1; etc
function SaveAxisSettings()
	string axisListStr, axisStr, tempStr
	variable numAxes, i
	
	string/G gAxisSettingsStr=""
	axisListStr = AxisList("") //list all axes on top graph
	numAxes = ItemsInList(axisListStr, ";")
	
	for (i=0; i<numAxes; i+=1)
		axisStr = StringFromList(i, axisListStr, ";")	//get one of the axis names
		gAxisSettingsStr +=axisStr						// add it to the list
		GetAxis/Q $axisStr								// get the current min and max of that axis
		sprintf tempStr, ";%0.5f;%0.5f;", V_Min, V_Max	//for now save 5 decimal places. Put this into a temporary string. Use sprintf to not lose digits like num2str would do.
		gAxisSettingsStr += tempStr						// append this string to the global string containing settings.
	endfor
	GetMarquee/K
end function

function RestoreAxisSettings()
	string axisListStr, axisStr
	variable numAxes, i, axisMin, axisMax, offsetinList

	if (exists("gAxisSettingsStr")!=2) //Test the global to see if it exists
		Abort "Global string 'gAxisSettingsStr' doesn't exist. Make sure to first run 'saveAxisSettings()' before running this function."
	else
		SVAR gAxisSettingsStr	//created in SaveAxisSettingsNEW
	endif
	axisListStr = AxisList("") //list all axes on top graph. (Could be different than list when SaveAxisSettingsNEW was run.)
	numAxes = ItemsInList(axisListStr, ";")
	
	for (i=0; i<numAxes; i+=1)	// loop thru all the current axes
		axisStr = StringFromList(i, axisListStr, ";")	//get an axis from current list of axes
		
// find axis named axisStr in saved settings gAxisSettingsStr - axisname0;min0;max0;axisname1;min1;max1; etc
		offsetinList = WhichListItem(axisStr, gAxisSettingsStr, ";", 0) // Find offset to that axis in the global string with the 'restore' settings.
		if (offsetinList !=- 1)	//if axis IS in string gAxisSettingsStr
			axisMin = str2num(StringFromList(offsetinList+1, gAxisSettingsStr, ";"))	//min is the next item in list after axis name
			axisMax = str2num(StringFromList(offsetinList+2, gAxisSettingsStr, ";"))	//max is 2 items away from axis name
			SetAxis $axisStr axisMin, axisMax											// set the axis with saved settings
		else	 // axis wasn't found in saved settings string
			Print "Axis "+axisStr+" not found in saved axis list. Autoscaling."
			SetAxis/A $axisStr
		endif
	endfor
	GetMarquee/K
end function

//Entering INF for the slope results in a vertical line with the intercept int. on the X axis.
macro MakeLineOnGraph(whichWave, XlineWaveName, YlineWaveName, slope, int)
	string whichWave, XlineWaveName, YlineWaveName
	variable slope, int
	prompt whichWave, "Select wave to draw line through:", popup,SortList(WaveList("*", ";", "WIN:"), ";", 16)
	prompt XlineWaveName, "Enter name for new X wave creating the line:"
	prompt YlineWaveName, "Enter name for new Y wave creating the line:"
	prompt slope, "Enter slope of line:"
	prompt int, "Enter Y intercept of line:"
	
	variable axisMin, axisMax
	string wInfo, whichSide, whichAxis, Xwave, Yaxis, wList
	
	if ((strlen(XlineWaveName)==0) + (strlen(YlineWaveName)==0))
		Abort "No wave name entered."
	endif
	
	wInfo = TraceInfo("", whichWave, 0)
	Yaxis = StringByKey("YAXIS", wInfo)	//pull out name of Y axis the selected wave is plotted on
	Xwave = StringByKey("XWAVE", wInfo)	//pull out name of X wave that selected wave is plotted against
	whichSide = StringByKey("AXISFLAGS", wInfo)	//pull out name of X axis that selected wave is plotted on

	GetAxis/Q bottom	//get the range of values for x axis in current graph
	axisMin = V_min
	axisMax = V_max

//	if ((int < axisMin) + (int > axisMax))
//		doAlert 0, "Line is out of the range of the selected wave."
//	endif
	XlineWaveName = CleanupName(XlineWaveName, 1 )
	YlineWaveName = CleanupName(YlineWaveName, 1 )

	Make/D/O/N=2 $XlineWaveName, $YlineWaveName

	if(numtype(slope)==1)	//if slope is inf
		$XlineWaveName =  {int, int}	//put y intercept into x wave
		if (cmpstr(whichSide, "/R")==0)	//Graph new line on right side
			getAxis/Q right					//get scale of axis that trace is on
		else
			getAxis left					//get scale of axis that trace is on
		endif
		$(YlineWaveName) = {V_min, V_max}	//set y wave to extent of left or right axis
	else
		$XlineWaveName =  {axisMin, axisMax}	//set x line to width of current graph x axis
		$(YlineWaveName) = {int, slope *( $(XlineWaveName)[1]-$(XlineWaveName)[0]) + int}	//calculate new line using x coordinates above
	endif
	
	wList = Wavelist("*", ";", "WIN:")
	if (FindListItem(YlineWaveName, wList, ";", 0)!=-1)	// See if the selected wave is on the graph already (will overwrite)
		RemoveFromGraph $YlineWaveName	//Remove it since user may have chosen a different wave to plot the line through but used existing line wave.
	endif

	if (cmpstr(whichSide, "/R")==0)	//Graph new line on right side
		AppendtoGraph/R=$YAxis $YlineWaveName vs $XlineWaveName
	else		//graph line on left side
		AppendtoGraph/L=$YAxis $YlineWaveName vs $XlineWaveName
	endif
end macro


menu "GraphMarquee"	//, dynamic	(commented out since it can cause delay with large graphs)
	"-"
	"Refresh axis list", BuildMenu "GraphMarquee"	//Build the Expand Sp menu on demand
	Submenu "Expand Sp"

		marqueeItem(0),/Q, ExpandSpecial(0)	//In Igor 6 when this returns an empty string, it doesn't show up in marquee menu.
		marqueeItem(1),/Q, ExpandSpecial(1)
		marqueeItem(2),/Q, ExpandSpecial(2)
		marqueeItem(3),/Q, ExpandSpecial(3)
		marqueeItem(4),/Q, ExpandSpecial(4)
		marqueeItem(5),/Q, ExpandSpecial(5)
		marqueeItem(6),/Q, ExpandSpecial(6)
		marqueeItem(7),/Q, ExpandSpecial(7)
		// as many as you think you'll ever have
	End
end menu

function/S marqueeItem(index)
	Variable index

	String topGraph= WinName(0,1)
	String menuStr="\\M0:(:_no more axes_"	// \\M0 to make disabling work on Mac and Win

	if( strlen(topGraph) )	//sometimes, command below produces list in different order.
		String aList= HVAxisList(topGraph,0)	// only left,right, etc. Only vertical axes.
		if (FindListItem("right", aList, ";", 0)!=-1)	//if this is in the list, pull it out of the list and make it 1st
			aList = "right;"+RemoveFromList("right", aList, ";")
		endif
		if (FindListItem("left", aList, ";", 0)!=-1)	//if this is in the list, pull it out of the list and make it 1st
			aList = "left;"+RemoveFromList("left", aList, ";")
		endif
		
		String thisAxis= StringFromList(index,aList)
		if( strlen(thisAxis) )
			menuStr=thisAxis
		endif
	endif
	return menuStr
End

Function ExpandSpecial(index)
	Variable index

	String topGraph= WinName(0,1)

	if( strlen(topGraph) )
		String aList= HVAxisList(topGraph,0)	// only left,right, etc, not bottom
		if (FindListItem("right", aList, ";", 0)!=-1)	//if this is in the list, pull it out of the list and make it 1st
			aList = "right;"+RemoveFromList("right", aList, ";")
		endif
		if (FindListItem("left", aList, ";", 0)!=-1)	//if this is in the list, pull it out of the list and make it 1st
			aList = "left;"+RemoveFromList("left", aList, ";")
		endif
		String thisAxis= StringFromList(index,aList)
		if( strlen(thisAxis) )
			GetMarquee/K $thisAxis
			SetAxis/W=$topGraph $thisAxis,  V_bottom, V_top
		endif
	endif
End

// Contributed by Harald Stark. Is NaN aware (which the built-in function StatsMedian is not).
Function calc_Median(w) // Returns median value of wave w.
	Wave w
	
	Variable result
	Duplicate/O w, tempMedianWave // Make a clone of wave
	Sort tempMedianWave, tempMedianWave // Sort clone
	SetScale/P x 0,1,tempMedianWave
	WaveStats/Q tempMedianWave
	result = tempMedianWave((V_npnts-1)/2)	//interpolates between values since brackets are round, not square
	KillWaves tempMedianWave
	return result
End


//Contributed by Harald Stark. Use this to trim a flight track which is overlayed onto an image. Sets points to NaN which
// are outside the range of the image. Removes the original flight track, and replaces it with the cropped one.
Function cropFlighttrack2Image()
	String im = ImageNameList("",";")

	if (itemsinList(im)>1)
		im = StringFromList(0,im)
		print "Warning! Found more than one image, will use:", im
	else
		im = StringFromList(0,im)
		print "using", im
	endif
	im = ReplaceString("'",im,"") // the first quotes contain the single quote to remove possible quotes from liberal names
	
	Wave w = $(im)
	
	Variable Latmin = min(DimOffset(w, 1),DimOffset(w, 1) + DimDelta(w,1)*DimSize(w,1))
	Variable Latmax = max(DimOffset(w, 1),DimOffset(w, 1) + DimDelta(w,1)*DimSize(w,1)) // sometimes scaling is negative
	Variable Lonmin = min(DimOffset(w, 0),DimOffset(w, 0) + DimDelta(w,0)*DimSize(w,0))
	Variable Lonmax = max(DimOffset(w, 0),DimOffset(w, 0) + DimDelta(w,0)*DimSize(w,0))
	String wlist = TraceNameList("",";",1)
	Variable i=0
	Variable stop = ItemsinList(wlist)
	String latname, lonname, latcropname, loncropname
	
	Do
		latname = NameOfWave(TraceNameToWaveRef("",StringFromList(i,wlist))) // This awkward looking construction avoids an error message when a wave is plotted twice in the graph
		lonname = XWaveName("", latname)
		if (strlen(lonname)>0)
			latcropname = Cleanupname(latname+"_cr",0)
			loncropname = Cleanupname(lonname+"_cr",0)
			Duplicate /O $latname, $latcropname
			Duplicate /O $lonname, $loncropname
			print "working on",latname,"vs.",lonname,", generated", latcropname,"and",loncropname
			Wave latcr = $latcropname
			Wave loncr = $loncropname
			RemoveFromGraph $latname
			AppendToGraph latcr vs loncr
			latcr = latcr>Latmin && latcr<Latmax ? latcr : NaN
			loncr = loncr>Lonmin && loncr<Lonmax ? loncr : NaN
		else
			print "could not work on", latname,", no x-wave found"
		endif 
		i+=1
	While (i<stop)
	legend
End Function

//------------------------------------------------------------------
// Macro makeFzSizeLegend() will create an f(z) size legend with as many symbols as there are sizes set in the f(z) settings. For example marker
// sizes 1-10 will result in a legend with 10 symbols. Marker sizes 5-10 will result in a legend with 6 symbols.
// Fractional symbol sizes aren't displayed on Mac screen, but they do print. PC's do seem to show the fractional size on screen.
// For even increments in legend, choose even bins for size. For example data 0-6000 and symbol size 5-10 will give 6 even bins:
// 0-1000, 1000-2000, etc up to 5000-6000.

// Numbers in the legend larger than 100000 and less than 0.01 will show up in scientific notation.

// At the prompt when macro starts, choose the wave that is displayed on the graph, not the f(z) wave.
// Works with liberal wave names.

// If you re-run the macro, the previous legend and waves plotted on the graph to make that legend are removed.
//Important: Needs to be run from datafolder containing f(z) wave (traceinfo doesn't give full datafolder if in a different data folder.)

macro makeFzSizeLegend(fzwave)
	string fzwave
	prompt fzwave, "Choose trace that uses f(z) wave from top graph: ", popup, sortlist(TraceNameList("", ";", 1), ";")

	makeSizeLegend(fzwave)
end

function makeSizeLegend(fzwave)
	string fzwave

	variable pos1, pos2, fz_min, fz_max, zMin, zMax, markerMin, markerMax, k, delta, markerNum
	variable numSymbols, legendMin, legendMax
	string fz_waveStr, fz_info, fz_infoSubStr,  s1, s2, legendPrecision, legendWaveStr, legendLabelStr, legendStr=""
		
	removeFzLegendWaves()	//clean up from previous f(z) legend if it is present.
	fz_info = traceInfo("", fzwave, 0)		// a long string containing all kinds of info including f(z) size info
	pos1 = strsearch(fz_info, "zmrkSize", 0)	//find "zmrkSize" starting at 0
	pos2 = strsearch(fz_info, ";", pos1+1)		// find first ; after zmrkSize
	fz_infoSubStr = fz_info[pos1, pos2]		//a substring starting with relevant f(z) info
	
	pos1=-1	//reset these for next search
	pos1 = strsearch(fz_infoSubStr, ",", 0)	//find "," in substring starting at 0. This should be the , following the wavename
	if (pos1 != -1)
		fz_waveStr=fz_infoSubStr[13, pos1-1]	//pull out wavename. Could be in ' ' if it's a liberal name
	else
		abort "Couldn't find comma denoting end of wave name."
	endif
	

//fz_infoSubStr should look like this now:
//zmrkSize(x)={GpsAlt,44,77,1,10}; 
//or
// zmrkSize(x)={'my liberal name',44,77,1,10}
//or 
//zmrkSize(x)={GpsAlt,*,*,1,10}; if auto size (one or both of the first two parameters can be *)
// use [A-Za-z] to read a string stopping when there's a character not within those brackets. (%s reads until white space which doesn't exist in this string)
// use [0-9.*e+-] to read a number or a * and possibly containing a decimal place and possibly scientific notation.

	sscanf fz_infoSubStr, "zmrkSize(x)={"+fz_waveStr+",%[0-9.*e+-],%[0-9.*e+-],%f,%f};", s1, s2, markerMin, markerMax

	fz_waveStr =ReplaceString("'", fz_waveStr, "" )	//remove ' if it is around a liberal name. (Even if it is, WaveStats doesn't want the ').
	zMin=NaN		//init. This is the range of the f(z) data wave.
	zMax=NaN

	Wavestats/Q $fz_waveStr	//Use wave min and max if not specified in zMin and zMax
	if (cmpstr(s1, "*")!=0)	//if there wasn't a * character, convert to a number
		sscanf s1, "%f", zMin	//convert s1 into number, allowing for long number that wouldn't translate properly with str2num
	else
		zMin = V_min
	endif
	if (cmpstr(s2, "*")!=0)
		sscanf s2, "%f", zMax	//convert s1 into number, allowing for long number that wouldn't translate properly with str2num
	else
		zMax = V_max
	endif
	
	pos1 = strsearch(fz_info, "marker(x)", 0)	//search later in info string for marker number that is plotted in graph.
	fz_infoSubStr = fz_info[pos1, pos1+20]	//get the substring that contains the symbol # used in graph
	sscanf fz_infoSubStr, "marker(x)=%d", markerNum
	
	numSymbols = 1+abs(markerMax-markerMin)	//could go in reverse order. +1 since we want both ends to show up on legend.
	legendMin = min(zMin, zMax)
	legendMax = max(zMin,zMax)
	delta = (legendMax - legendMin) / (numSymbols)

//figure out how many decimal places are needed to display the numbers in the legend.
	if ((delta >100000) + (delta < 0.01))	//use scientific notation outside of this range.
		legendPrecision = "%0.2e"
	elseif (delta < 0.1)
		legendPrecision = "%0.3f"
	elseif (delta < 1)
		legendPrecision = "%0.2f"
	elseif (delta < 10)
		legendPrecision = "%0.1f"
	else
		legendPrecision="%0.0f"
	endif

	
	String savedDF= GetDataFolder(1)	// Remember CDF in a string.
	NewDataFolder/O/S root:fzLegendStuff	//make new datafolder to contain waves we need for legend.
	Killwaves/A/Z							//kill all existing waves in this data folder
	for (k=0; k < numSymbols; k+=1)		//loop over number of bins for legend.
		legendWaveStr = "legendWave"+num2str(k)	//Make waves called legendWave0, legendWave1, etc
		Make/o/d/n=1 $legendWaveStr=NaN				//create waves and set to NaN

		sprintf legendLabelStr, legendPrecision, legendMin + (k*delta) + delta/2	//safely convert num 2 string without dropping digits in large numbers.
		legendStr +="\s("+legendWaveStr+") "+legendLabelStr+" \r"	//build the legend string
		AppendToGraph/L=hiddenLeft/B=hiddenBottom $legendWaveStr	//append these waves for use in the legend.
		ModifyGraph marker($legendWaveStr)=markerNum		// set the marker for the legend to same marker used in f(z)
		if (markerMax > markerMin)
			ModifyGraph mode($legendWaveStr)=3, msize($legendWaveStr)=k+markerMin	// increasing order is used for symbol size.
		else
			ModifyGraph mode($legendWaveStr)=3, msize($legendWaveStr) = markerMin-k		//reverse order is used.
		endif
	endfor

	legendStr=RemoveEnding(legendStr)	//chop off last \r
	ModifyGraph noLabel(hiddenLeft)=2,noLabel(hiddenBottom)=2,axThick(hiddenLeft)=0	//hide the axes used to plot the legendwaves
	ModifyGraph axThick(hiddenBottom)=0
	Legend/C/N=fz_LegendText/J/M legendStr	//put the size legend on the graph, sized by symbol size in graph.

	SetDataFolder savedDF	// Restore CDF from the string value.
end

function removeFzLegendWaves()
	variable n, num, DFexists
	string wlist, wStr
	String savedDF= GetDataFolder(1)	// Remember CDF in a string.

	DFexists = DataFolderExists("root:fzLegendStuff")
	if (DFexists)	//else do nothing since none of this exists
		SetDataFolder root:fzLegendStuff	//make new datafolder to contain waves we need for legend.
		
		Legend/K/N=fz_LegendText	//kill legend if it exists.
		wlist = WaveList("legendWave*", ";", "WIN:" )
		num = ItemsInList(wlist, ";")
		for (n=0; n<num; n+=1)
			wStr = StringFromList(n, wlist, ";")
			RemoveFromGraph $wStr
		endfor
		SetDataFolder savedDF	// Restore CDF from the string value.
	endif
end

// Code contributed by Harald Stark. Calculates average vertical profile using either average or median.
// Outputs waves with "_ap" appended on the end for the averaged waves. Also outputs N, sdev waves to be used as error bars.
Macro calc_avg_alt_profile(datwav, altwav, startalt, stopalt, altbinsize, startpt, stoppt, levelonly, method)
	String datwav
	String altwav = "GPSAlt"
	Variable startalt = 0
	Variable stopalt = 5000
	Variable altbinsize = 100
	Variable startpt = -1
	Variable stoppt = -1
	String levelonly = "No"
	String method = "Average"
	Prompt datwav,"Select data wave", popup SortList(WaveList("*",";",""), ";", 16)
	Prompt altwav,"Select altitude wave", popup SortList(WaveList("*",";",""), ";", 16)
	Prompt startpt, "Start pt. num (-1 for beginning)"
	Prompt levelonly, "Use leveled flight parts only?",popup "Yes;No"
	Prompt method, "Method",popup "Average;Median"
	
	if (numpnts($datwav)!=numpnts($altwav))
		print datwave,"and",altwave,"don't have the same number of points!"
	else
		Do_calc_avg_alt_profile(datwav, altwav, startalt, stopalt, altbinsize, startpt, stoppt, levelonly, method)
		String df = GetDataFolder(1)
		print "data folder:", df
	endif

End Macro

// Does the work for the above macro
Function Do_calc_avg_alt_profile(dwstr, awstr, startalt, stopalt, altbinsize, startpt, stoppt, levelonly, method)
	String dwstr, awstr, levelonly, method
	Variable startalt, stopalt, altbinsize, startpt, stoppt
	
	if (startpt>0) // start point selected
		if (stoppt<0) // stop point set to end of wave
			stoppt = numpnts($dwstr) -1
		endif
		duplicate/O/R=[startpt,stoppt] $dwstr, dwtempwav
		duplicate/O/R=[startpt,stoppt] $awstr, awtempwav
	else // no start or stop point selected, duplicate whole wave
		duplicate/O $dwstr, dwtempwav
		duplicate/O $awstr, awtempwav
	endif
	if (cmpstr(levelonly,"Yes")==0)
		if (exists("roll")&&exists("pitch"))
			wave dwtemp = $("dwtempwav")
			wave rol = $("roll")
			wave pit = $("pitch")
			dwtemp = (abs(rol)<3&&abs(pit)<3) ? dwtemp : NaN
		else
			abort "can't find roll and pitch waves"
		endif
	endif
	
	Sort awtempwav dwtempwav,awtempwav
	String alt_prof = cleanupname(dwstr + "_ap",0)
	if (cmpstr(method,"Median")==0)
		alt_prof = cleanupname(dwstr + "_apm",0)
	endif
	String alt_prof_sdev = cleanupname(dwstr + "_ap_sdev",0)
	String alt_prof_N = cleanupname(dwstr + "_ap_N",0)
	if (cmpstr(alt_prof, alt_prof_sdev)==0)
		alt_prof = cleanupname(dwstr[0,27] + "_ap",0)
		alt_prof_sdev = cleanupname(dwstr[0,25] + "_ap_sd",0)
	endif
	// String prof_alt = cleanupname(dwstr + "_profalt",0)
	String prof_alt = cleanupname("alt_"+num2str(startalt)+"_"+num2str(stopalt)+"_"+num2str(altbinsize),0) 
	Variable numalts = ceil(stopalt - startalt) / altbinsize
	Make/N=(numalts)/O $alt_prof, $alt_prof_sdev, $prof_alt, $alt_prof_N
	SetScale/P x (startalt+altbinsize/2),(altbinsize),"", $alt_prof, $alt_prof_sdev, $prof_alt, $alt_prof_N
	if (startpt>0) // start point selected
		Print "Created profile from [",startpt,",",stoppt,"]:", alt_prof, ", profile altitudes:", prof_alt, ", standard deviations:", alt_prof_sdev
	else
		Print "Created profile:", alt_prof, ", profile altitudes:", prof_alt, ", standard deviations:", alt_prof_sdev
	endif
	Wave ap = $alt_prof
	Wave aps = $alt_prof_sdev
	Wave apn = $alt_prof_N
	Wave pa = $prof_alt
	ap = nan
	aps = nan
	apn = nan
	// calculate altitude wave
	// pa = startalt + altbinsize / 2 + p * altbinsize
	
	pa = nan
	Variable start = FindFirstGE(awtempwav, startalt, 0)
	Variable j= 0 //floor((awtempwav[0]-startalt)/altbinsize) // starting index
	Variable stop = numpnts(awtempwav)
	Variable binstop = startalt + altbinsize//altbinsize*ceil(awtempwav[0]/altbinsize) // set to top of lowest alt bin
	Variable currpt = start
	Variable medi
	Variable i = start
	Do
		i = FindFirstGE(awtempwav, binstop, i) // set index to start of next bin (or end of wave)
		
		Wavestats/R=[currpt,i-1]/Q dwtempwav // get avg data value for all points below top of bin
		if (cmpstr(method,"Median")==0)
			medi = Median_StartStop(currpt,i-1, dwtempwav)
			ap[j] = medi
		else
			ap[j] = V_avg
		endif
		//NVAR Va = V_avg
		//NVAR Vs = V_sdev
		
		aps[j] = V_sdev
		apn[j] = V_npnts
		Wavestats/R=[currpt,i-1]/Q awtempwav // get avg alt value for all points below top of bin
//		PRINT "J, currpt, I-1, V_AVG= ", J, currpt, I-1, V_AVG
		pa[j] = V_avg
		binstop += altbinsize
		currpt = i
		j+=1
//		print binstop, stopalt, stop
	While (binstop<=stopalt&&i<stop)
	Killwaves awtempwav, dwtempwav
End Function

//The built-in Igor function for median can't handle NaNs. This routine can.
Function Median_StartStop(start, stop, w) // Returns median value of wave w, points between start and stop, Nan's can be present
	Variable start, stop
	Wave w
	
	Variable result
	Duplicate/O/R=[start,stop] w, tempMedianWave // Make a clone of wave
	Sort tempMedianWave, tempMedianWave // Sort clone
	SetScale/P x 0,1,tempMedianWave
	WaveStats/Q tempMedianWave
	result = tempMedianWave((V_npnts-1)/2)
	// This following line was in the original function, it will give wrong results when Nan's are present in the wave
	//result = tempMedianWave((numpnts(tempMedianWave)-1)/2)
	KillWaves tempMedianWave
	return result
End

// Makes legend with name of wave and axis it is plotted on. Labels all axes too. Handy for very busy graph where you can't tell which axis
// is used for a trace. Running the function again turns off all the labels that were created by running it.
// Taken from the Igor mail list.
Function toggleaxisinfo() // for each axis -> show name as tag
	string al=AxisList(""), an, ai // axis list, axis name, axis info
	variable k,n=ItemsInList(al),gol,mid,vh
	variable ta=strlen(GetUserData("","","ax_tag")) // tag on axis

	GetMarquee/K
	for(k=0;k<n;k+=1)
		an=StringFromList(k,al)
		if (ta) // remove tags
			Tag/N=$(an)/K
			TextBox/K/N=ax_cap
			SetWindow kwTopWin UserData($"ax_tag")=""
		else // append tags
			ai=AxisInfo("",an)
			gol=NumberByKey("log(x)", AxisInfo("",an),"=")
			GetAxis/Q $(an) // V_min, V_max
			if (gol)
				mid=alog((log(V_min)+log(V_max))/2)
			else
				mid=(V_min+V_max)/2
			endif
			vh=Vertical(an)?90:0 // label orientation
			Tag/N=$(an)/C/F=0/O=(vh)/B=(65535,49151,49151)/X=0/Y=0/L=0 $(an), mid, an
			caption()
			SetWindow kwTopWin UserData($"ax_tag")="tag"
		endif
	endfor
end

Function Vertical(ax_name) // is the axis vertical ?		//Used by function toggleaxisinfo()
	String ax_name // 0 = horizontal, 1 = vertical, -1 = not present
	String side=AxisInfo("",ax_name), car="top;left;bottom;right"

	Return mod(WhichListItem(StringByKey("AXTYPE",side),car),2)
End // return -1 when axis is not present on top graph

Function caption() // "trace -> X, Y axis"	//Used by function toggleaxisinfo()
	string al=TraceNameList("",";",1), an, ti,tx
	variable k,n=ItemsInList(al)

	tx="trace -> \tX, Y axis"
	for(k=0;k<n;k+=1)
		an=StringFromList(k,al)
		ti=TraceInfo("",an,0)
		tx+="\r\s("+an+") "+an+" \t"
		tx+=StringByKey("XAXIS", ti)+", "
		tx+=StringByKey("YAXIS", ti)
	endfor
	TextBox/C/N=ax_cap/B=(65535,54611,49151)/A=MC/X=0/Y=0 tx
end

// Set values in a keywaveNaN to NaN at the transitions in keywave.
// Specify number of points to NaN after the transition: numAfterTransitionToNaN
// Specify whether to use NaN as the transition: useNaNsAsTransitions. Ignores NaNs if set to 0
function NanTransitions(keywaveStr, keywaveToNaNStr, numAfterTransitionToNaN, useNaNsAsTransitions)
	string keywaveStr, keywaveToNaNStr
	variable numAfterTransitionToNaN, useNaNsAsTransitions
	
	variable n, len
	
	wave keywave=$keywaveStr
	wave keywaveNaN=$keywaveToNaNStr
	len = numpnts(keywave)

//	keywaveNaN = keywave
	for(n=0; n < len; n+=1)
		if (useNaNsAsTransitions)
			if (numtype(keywave[n]) != numtype(keywave[n+1]))
				keywaveNaN[n+1, n+numAfterTransitionToNaN] = NaN
			endif
		else
			if ((numtype(keywave[n])==0) * (numtype(keywave[n+1])==0))
				if (keywave[n] != keywave[n+1])
					keywaveNaN[n+1, n+numAfterTransitionToNaN] = NaN
				endif
			endif	//else ignore the NaN if they aren't being used as a transition
		endif
	endfor
end

// User selects points in a scatter plot using the marquee.
// All those points are then copied into wave "tracewHighlighted" and then plotted in a time series.
// If the trace uses f(z) then there's an option to highlight the f(z) data instead of the Y data.
// Shows where points are in time series.
// KA 20110124
//--------- beginning of code for FindPointsInTimeSeries()
Function FindPointsInTimeSeries() : GraphMarquee
	string xwname, wlist, wstr, timewaveStr, yaxis, fzWave
	variable ind, len, newgraphFlag, bottomVal, topVal, leftVal, rightVal
	
	DoWindow/F ScatterplotToTimeSeriesPanel	//Bring panel to front
	if (!V_flag)	//panel doesn't exist
		Execute "ScatterplotToTimeSeriesPanel()"	//Make the panel
		DoAlert 0, "Please use panel to make settings before running FindPointsInTimeSeries(). Leave panel open after making selections."
		abort	//Do nothing if panel is not found.
	endif

	ControlInfo/W=ScatterplotToTimeSeriesPanel Tracepopup	//get menu selection for Y wave
	if (V_Flag==0)
		DoAlert 0, "Couldn't find popup menu in panel."
	endif
	if (cmpstr(S_Value, "none")==0)							//user didn't make selection in menu, and menu is still set on default value
		Abort "Please choose Y wave from menu."
	endif
	wstr = S_Value												// Y trace to hightlight

	yaxis = StringByKey("YAXIS", TraceInfo("",wstr, 0))	//get name of Y axis that selected trace is plotted on (ie. left, right)

	GetMarquee $yaxis, bottom	//get marquee info. Assumes data on graph is using x (bottom) and y (which ever side) axes
	if (V_Flag == 0)
		Print "There is no marquee"
		return 0
	endif
	bottomVal = V_bottom	//save the values from GetMarquee
	topVal = V_top
	leftVal = V_left
	rightVal = V_right
	
	fzWave=""				//init
	ControlInfo checkFZ		//See if the checkbox is checked
	if (V_value)
		fzWave = extractFZinfo(wstr)	// get name of the f(z) wave
		if (strlen(fzWave)>0)			// Found name 
			Print "f(z) wave selected is ", fzWave
			wave fz = $fzWave
		else
			abort "No f(z) data is plotted."	//Box is checked, but graph doesn't use f(z) trace.
		endif
	endif
	
	wave tracew = TraceNameToWaveRef("", wstr )	//the wave name of the trace selected
	xwname = XWaveName("", wstr)				//name of the xwave that first trace is plotted against
	if (exists(xwname)==1)						//Make sure that wave exists.
		wave xw = $xwname
	else
		abort "Couldn't find wave "+xwname+" in current data folder."
	endif

	ControlInfo/W=ScatterplotToTimeSeriesPanel Timewavepopup	//get menu selection for time wave
	if (V_Flag==0)
		DoAlert 0, "Couldn't find check box in panel."
	endif
	if (cmpstr(S_Value, "none")==0)								//user didn't make selection in menu, and menu is still set on default value
		Duplicate/O $wstr, TempXwave								//create a temporary wave to function as "calculated"
		TempXwave=p												// = point number
		timewaveStr = "TempXwave"								// set string to TempXwave so rest of function will work.
//		Abort "Please choose time wave from menu."					// Don't abort anymore if no selection for timewave
	else
		timewaveStr = S_Value										// name of timewave selected in menu
	endif

	Print "Finding points in waves ", wstr, "and x wave ", xwname+" . Timewave is ", timewaveStr
	len =  numpnts(tracew) 									//number of points in trace
	Make/O/D/N=(len) tracewHighlighted=NaN					//This wave will contain the points within the marquee 
	
	for (ind=0; ind<len; ind+=1)								// loop through all points in trace
		if ((tracew[ind] > bottomVal) && (tracew[ind] < topVal) && (xw[ind] > leftVal) && (xw[ind] < rightVal))	//within marquee
			if (strlen(fzWave))
				tracewHighlighted[ind] = fz[ind]	//found point in f(z) data within marquee. Copy it to tracewHighlighted wave.
			else
				tracewHighlighted[ind] = tracew[ind]	//found a point in Y data trace within marquee. Copy it to tracewHighlighted wave
			endif
		endif
	endfor
	
	ControlInfo/W=ScatterplotToTimeSeriesPanel newGraphCheck			//read checkbox for creating new graph
	if (V_Flag==0)
		DoAlert 0, "Couldn't find check box in panel"
	endif
	newgraphFlag = V_Value
	if (newgraphFlag)														//If checkbox is checked for making a new graph
		display tracewHighlighted vs $timewaveStr 							// display graph of highlight trace
		ModifyGraph mode=3, rgb(tracewHighlighted)=(65535,65535,0)	// yellow markers
		ModifyGraph marker(tracewHighlighted)=19,msize(tracewHighlighted)=6	//solid circle for marker, size 6
		if (strlen(fzWave))
			AppendToGraph fz vs $timewaveStr								// append the whole f(z) wave that was on the first graph
		else
			AppendToGraph tracew vs $timewaveStr								// append the whole Y trace that was on the first graph.
		endif
	endif
End

// TraceInfo returns either
//   "zColor(x)=0" for no f(z) selected
// or
//   "zColor(x)={CO2_ppmv,*,*,Rainbow,1}" for f(z) color-coding
// This function will return the f(z) wave name by extracting it from the string which TraceInfo returns
function/S extractFZinfo(yWave)
	string yWave
	
	string teststr, fzStr=""
	variable pos1, pos2
	
	teststr = StringByKey("RECREATION", TraceInfo ("", yWave, 0))	//get the string that follows "RECREATION" which has the f(z) info.
	
	if (cmpstr(teststr, "zColor(x)=0")!=0)	//if no f(z) is selected for this trace, then string = "zColor(x)=0"
		pos1 = strsearch(teststr, "{", 0)		//find position of "{"
		pos2 = strsearch(teststr, ",", pos1)		// find position of "," which follows f(z) name
		fzStr = teststr[pos1+1, pos2-1]		//extract name of f(z) wave.
	endif
	return fzStr
end

Window ScatterplotToTimeSeriesPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(395,50,737,366)
	ModifyPanel frameInset=2
	SetDrawLayer UserBack
	SetDrawEnv fsize= 9
	DrawText 56,225,"Makes a time series of the highlighted data"
	SetDrawEnv fsize= 10
	DrawText 26,54,"Top graph needs to be the one containing data to highlight."
	SetDrawEnv fsize= 10
	DrawText 26,67,"Select data with marquee."
	PopupMenu Timewavepopup,pos={42,143},size={206,20},title="Time wave for data in scatterplot:"
	PopupMenu Timewavepopup,mode=3,popvalue="none",value= #"WaveList(\"*\",\";\", \"\")"
	CheckBox newGraphCheck,pos={54,193},size={149,14},title="New graph of highlighted data"
	CheckBox newGraphCheck,value= 0
	GroupBox group0,pos={43,172},size={260,66},title="Output Options",frame=0
	GroupBox group1,pos={18,21},size={308,257},title="Settings for \"Scatterplot to time series\""
	GroupBox group1,fSize=10,fStyle=0
	PopupMenu Tracepopup,pos={43,83},size={145,20},title="Trace to highlight:"
	PopupMenu Tracepopup,mode=1,popvalue="none",value= #"TraceNameList(\"\",\";\", 1)"
	CheckBox checkFZ,pos={44,111},size={184,14},title="Highlight f(z) wave instead of Y wave."
	CheckBox checkFZ,value= 0
EndMacro
//---- end of code for FindPointsInTimeSeries()



// *************************************************
//  Convert 1D waves to 2D matrix functions
// *************************************************


// *************************************************
// from a table, creates a 2D matrix of all waves in the table that has the name wavePrefixblah
// the name of the newly created 2D wave is wavePrefixMx
// if no TabelNameStr is supplied, code will use top table
// example useage Create2Dfrom1DFromTable("Org_dMdlogDva_", 1)
Function Create2Dfrom1DFromTable(wavePrefix, makeTableFlag,[TableNameStr])
string wavePrefix, TableNameStr	
variable makeTableFlag

	variable idex, numrow, numcol, rowVal, colVal, prefixLen
	string myWaveList, myWaveStr, mytableInfo
	
	if (ParamIsDefault(TableNameStr ))
		TableNameStr=""
	endif
	
	prefixLen=strlen(wavePrefix)
	mytableInfo = TableInfo(TableNameStr, -2)
	numrow = str2num(StringByKey("Rows",mytableInfo))
	numcol = str2num(StringByKey("Columns",mytableInfo))
	
	Make/o/d/n=(numrow, numcol) $wavePrefix+"Mx"/wave=my2Dmatrix
	
	rowVal=0
	colVal=0
	for(idex=0;idex<numcol;idex+=1)
		mytableInfo = TableInfo(TableNameStr, idex)
		myWaveStr =  StringByKey("COLUMNNAME",mytableInfo)
		if (stringmatch(myWaveStr[0,prefixLen-1], wavePrefix) )
			wave wav = $StringByKey("Wave",mytableInfo)
			if (rowVal==0)
				rowVal=numpnts(wav)
			endif
			if (numpnts(wav)!= rowVal || dimsize(wav,2)>0)
				abort "something was wrong with wave "+myWaveStr
			endif
			
			ImageTransform/G=(colVal)/D=wav putCol my2Dmatrix
			colVal+=1
		endif
	endfor

	redimension/n=(RowVal,ColVal) my2Dmatrix
		
	if (makeTableFlag)
		edit my2DMatrix
	endif

End


// *************************************************
// use this if the wave names can be sorted alphabetically (i.e. name_01, name_02) and this order is the order you want
// example usage  Create2Dfrom1DFromDataFolder("Org_dMdlogDva_", 1)
Function Create2Dfrom1DFromDataFolder(wavePrefix,makeTableFlag)
string wavePrefix
variable makeTableFlag

	variable idex, numrow, numcol, rowVal, colVal, prefixLen
	string myWaveList, myWaveStr, mytableInfo
	
	prefixLen=strlen(wavePrefix)
	myWaveList = WaveList(wavePrefix+"*", ";", "DIMS:1" )

	myWaveList = SortListAlphabetically(myWaveList)
	wave wav = $stringFromList(0, myWaveList)
	numrow = numpnts(wav)
	numcol = itemsInList(myWaveList)
	
	Make/o/d/n=(numrow, numcol) $wavePrefix+"Mx"/wave=my2Dmatrix
	
	for(idex=0;idex<numcol;idex+=1)
		wave wav = $stringFromList(idex, myWaveList)
		ImageTransform/G=(colVal)/D=wav putCol my2Dmatrix
	endfor

	if (makeTableFlag)
		edit my2DMatrix
	endif
End


// Puts a colon on the end of a string if it doesn't have one or does the reverse if you set none=1
Function/S AddOrRemoveColon(Str [,none])
	String Str
	Variable none
	
	if (ParamIsDefault(none))		// Put colon on end
		if (!stringmatch(str, "*:"))
			str += ":"
		endif
	else						// Take colon off end
		if (stringmatch(str, "*:"))
			String endstr=""
			Variable n
			for (n=0; n<ItemsInList(str, ":"); n+=1)
				endstr  += StringFromList(n, str, ":")
				if (n<(ItemsInList(str, ":") - 1))
					endstr += ":"
				endif
			endfor
			str = endstr
		endif
	endif

	return str
End
