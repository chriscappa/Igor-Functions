#pragma rtGlobals=1		// Use modern global access method.

//*************************************************************************************************************************************************
Function TopOfSizeDist_toolkit_Procs() 
	//The only purpose of this function is to be the reference point for opening this procedure window from the menu bar.
	//=> leave this function at the top of the procedure window
End
//*************************************************************************************************************************************************

//SMPS Loading and Processing Functions
//	By Chris Cappa (UC Davis) building upon Tim Onasch's (Aerodyne) SMPS code
//
//****Purpose:  Load and process SMPS data obtained with TSI SMPS software and "Exported to File" as *.txt file 
//	   with data type = dW/dlogDp, Number; Delimiter = Tab; Orientation = ROWS.  The loading program either loads
//	   a single file or all *.txt files (unless explicitly included in exclusion list) into root:smps: data folder.  Can add as
//	   many single or multiple files as desired and will PROCESS all of the loaded data into one set of data under
//	   root:smps: data folder.  Can now post-process in any data folder (e.g. root:smps2:) CANNOT as yet deal with 
//	   multiple files with varying diameter bases!
//*****Purpose2: Load and Process APS data obtained with TSI APS software and "Exported to File" as *.txt file with either
//		data type = dW/dlogDp or dW, Number; Delimiter = Tab; Orientation = ROWS. 
//
//ToDoList working on:
//0.)  Allow program to work with dissimilar diameter wave data
//1.)  Use liberal names for file names of loaded SMPS data (smps.txt)
//2.)  NOTE: TSI program does charge correction incorrectly!!!!!!!!
//3.)  Sample time is actually just tdiff - which causes problems when user waits between scans.  TSI does not export scan time info for some stupid reason.
//
//
//****ToDoList
//	1.  Deal with files that have differing diameter scales.  Either via a simple check and warning message, or a true fix.
//	2.  Be smart about choosing files in folder by loading and reading first line or so.
//
//****Function INDEX --> NEEDS TO BE UPDATED
//	load_single_smps_file() - allows user to identify single *.txt smps file to load
//	function load_smps_directory() - allows user to identify folder with single or multiple *.txt smps files to load
//	function load_smps_row(file_list) - loads all smps *.txt files listed in file_list into subfolders in data folder root:smps:
//	function process_smps_data() - processes all data in all subfolders in root:smps: and concatenates results into root:smps:


//****Data format for SMPS data - can be loaded with/without header or exact header
//	Sample File	C:\DOCUME~2\Aerodyne\Projects\SP2\BOSTON~1\SMPS\200710~1.S80
//	Classifier Model	3080
//	DMA Model	3081
//	DMA Inner Radius(cm)	0.00937
//	DMA Outer Radius(cm)	0.01961
//	DMA Characteristic Length(cm)	0.44369
//	CPC Model	3010
//	Gas Viscosity (kg/(m*s))	1.8203e-005
//	Mean Free Path (m)	6.65e-008
//	Channels/Decade	64
//	Multiple Charge Correction	TRUE
//	Nanoparticle Aggregate Mobility Analysis	FALSE
//	Diffusion Correction	TRUE
//	Gas Density	0.0012
//	Units	dw/dlogDp
//	Weight	Number
//	Sample #	Date	Start Time	Diameter Midpoint	7.64	7.91	8.20	8.51	8.82	9.14	9.47	9.82	 10.2	... truncated here due to line length issues
//	1	10/11/07	13:48:31		1.22837e+006	939869	1.03757e+006	897689	888182	947596	... truncated here due to line length issues
//	2	10/11/07	13:50:35		3.48917e+006	2.47102e+006	2.00932e+006	1.56242e+006	... truncated here due to line length issues

menu  "SMPS and APS"
	
	"Show Size Dist procs"
	
	Submenu "Data Loading"
		"Load SMPS Files", load_smps_files()
		"Load APS Files", load_aps_files()
		"Load CPC Files", load_cpc_files()
		"Load SEMS Files", Load_SEMS_files()
	end
	Submenu "Averaging"
		"AveragSizeDistOnlyXSecs"
		"AveragSizeDistUsingStartStop"
	end
	Submenu "TD"
		"SplitSMPS_byTDstate", SplitSMPS_byTDstate()	
		"Average between cycles", InterpSMPSforTD()
		"SplitSEMS_byTDstate", SplitSEMS_byTDstate()
	end
	Submenu "Data Merging"
		"MergeSMPS_tseries", MergeSMPS_tseries("dNdlogDp")
		"MergeAPS_tseries", MergeAPS_tseries("dN")
		"MergeCPC_tseries", MergeCPC_tseries()
		"MergeSEMS_tseries", MergeSEMS_tseries("dNdlogDp",0)
		"MergeSMPSwithAPS", MergeSMPSwithAPS(600,2400)
	end
	"ConvertDpaToDpm - Constant Density", ConvertDpaToDpm()
	"ConvertDpaToDpm - Density Wave", ConvertDpaToDpm_DensityWave()
	"dNdlogDp to Surface and Volume", Convert_dNdlogDp(1)
	"Integrate Size Distributions",  IntegrateSizeDist()
end

Function ShowSizeDistprocs()
	DisplayProcedure "TopOfSizeDist_toolkit_Procs"
End

//===================Names============
// sdm = size_dist_matrix
// dt = datetime_wave

//===================User defined variables=======================
//Path to local data source
//strconstant local_path_to_data = "C:\documents\Aerodyne\AMS field\Mexico City 2006\PNNL CTOF\SMPS"
//constant density=1.0
//===================User defined variables=======================

//************************************************************************************
function load_smps_files()
	// to load a single SMPS file
	string extension = "txt"
	Variable refNum
	String message = "Select one or more SMPS files"
	String file_list
	String fileFilters = "Data Files (*."+extension+"):."+extension+";"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	file_list = S_fileName
	if (strlen(file_list) == 0)
		print "==================================="
		print "No Data Files Chosen for Loading.  Aborting."
		print "==================================="
		Close/A
		setdatafolder root:
		abort
	endif
	Close/A
	file_list = replacestring("\r",file_list,";")
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:smps	
	load_smps_row(file_list)
	setdatafolder root:smps
end

//************************************************************************************
function load_smps_directory()
	// to load all SMPS files in a directory
	newpath/o/q path_to_data	//, local_path_to_data
	string file_list = IndexedFile(path_to_data, -1, ".txt")
	if (cmpstr(file_list, "")==0)
		print "==================================="
		print "No SMPS *.txt data found in directory.  Aborting."
		print "==================================="
		setdatafolder root:
		abort
	endif
	//Enter *.txt filenames to remove from list (i.e. SKIP PROCESSING) using format in line below
//	file_list = RemoveFromList("LongDMA_calnex_SanDiego_11May2010_spl_400.txt", file_list)
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:smps	
//	return file_list	
	load_smps_row(file_list)
	setdatafolder root:smps
end

function load_smps_row(file_list)
	string file_list 	//= initialize_smps()
	string sdf = getdatafolder(1)
	string file_full_path, file_name, file_name_str
	variable nf = itemsinlist(file_list)
	variable i_nf, nc, i_header, i
	for (i_nf=0;i_nf<nf;i_nf+=1)
		file_full_path = stringfromlist(i_nf,file_list)
		file_name = ParseFilePath(0, file_full_path, ":", 1, 0)
		file_name_str = ParseFilePath(3, file_full_path, ":", 0, 0)

		 if (strlen(file_name)>31)
//		 	file_name = file_name[0,30]
		 endif
		 
	//Test if the filename starts with number or not
		 if (numtype( str2num(file_name[0]) ) ==0)
		 	file_name_str = ("root:smps:x"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 else
		 	file_name_str = ("root:smps:"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 endif
		newdatafolder/o/s $file_name_str
		LoadWave/J/K=2/V={""," $",0,0}/L={0,0,0,0,1}/n=smps file_full_path
		wave/t smps0
		variable nr=numpnts(smps0)
	//Header INFO: Test if file contains header info, if so count header lines for appropriate wave, and Load Header Parameters
		if (numtype(str2num(smps0[0]))!=0)
			i_header=0
			for (i=0;i<nr;i+=1)
				if (numtype(str2num(smps0[i]))!=0)
					i_header+=1
				else
					break
				endif
			endfor
			//hdr = header info
			make/o/n=(i_header-1,2)/t  $("hdr_"+stringfromlist(0,S_fileName,"."))
			wave/t smps_header_info = $("hdr_"+stringfromlist(0,S_fileName,"."))
			smps_header_info[][0] = stringfromlist(0,smps0[p],"\t")
			smps_header_info[][1] = stringfromlist(1,smps0[p],"\t")
		endif			
		
	//Load Size Distribution Waves
		string d,t
		variable np=0,d1,d2,d3,t1,t2,t3
		//sn = sample number
		make/o/d/n=0 $("sn_"+stringfromlist(0,S_fileName,"."))
		wave sample_no = $("sn_"+stringfromlist(0,S_fileName,"."))
		make/o/d/n=0 $("dt_"+stringfromlist(0,S_fileName,"."))
		wave dt = $("dt_"+stringfromlist(0,S_fileName,"."))
		make/o/d/n=0 $("tdiff_"+stringfromlist(0,S_fileName,"."))
		wave tdiff = $("tdiff_"+stringfromlist(0,S_fileName,"."))
		for (i=0;i<nr;i+=1)
			if (numtype(str2num(smps0[i]))==0)
				insertpoints np,1,dt,sample_no,tdiff
				sample_no[np]=str2num(stringfromlist(0,smps0[i],"\t"))
				d=stringfromlist(1,smps0[i],"\t")
				d1=str2num(stringfromlist(0,d,"/"))
				d2=str2num(stringfromlist(1,d,"/"))
				d3=str2num(stringfromlist(2,d,"/"))
				if(d3<2000)
					d3+=2000	// in case year loads as YY instead of as YYYY
				endif
				t=stringfromlist(2,smps0[i],"\t")
				t1=str2num(stringfromlist(0,t,":"))
				t2=str2num(stringfromlist(1,t,":"))
				t3=str2num(stringfromlist(2,t,":"))
				dt[np]=date2secs(d3,d1,d2)+t1*3600+t2*60+t3
				tdiff[np] = dt[np] - dt[(np-1)]
				if (np==1)
					tdiff[0] = tdiff[np]
				endif
				np+=1
			endif
		endfor
		make/o/d/n=0 $("dia_"+stringfromlist(0,S_fileName,"."))
		wave diameter = $("dia_"+stringfromlist(0,S_fileName,"."))
		variable nd,j,k=0,l,m,start_sdm=0
		for (i=0;i<nr;i+=1)
			if (stringmatch(stringfromlist(0,smps0[i],"\t"),"Sample #"))
				nd=itemsinlist(smps0[i],"\t")
				for (j=0;j<nd;j+=1)
					if (numtype(str2num(stringfromlist(j,smps0[i],"\t")))==0)
						insertpoints k,1,diameter
						diameter[k]=str2num(stringfromlist(j,smps0[i],"\t"))
						k+=1
						if(start_sdm==0)
							start_sdm = j
						endif
					endif
				endfor
				break
			endif
		endfor
		//sdm = size distribution matrix
		make/o/d/n=(0,k) $("sdm_"+stringfromlist(0,S_fileName,"."))
		wave sdm = $("sdm_"+stringfromlist(0,S_fileName,"."))
		m=0
		l=0
		for (i=0;i<nr;i+=1)	
			if (numtype(str2num(smps0[i]))==0)
				nd=itemsinlist(smps0[i],"\t")
				insertpoints l,1,sdm
				for (j=start_sdm;j<k+start_sdm;j+=1)
					sdm[l][m]=str2num(stringfromlist(j,smps0[i],"\t"))
					m+=1
				endfor
				l+=1
				m=0
			endif
		endfor
		// turn NaN values into 0 values
		sdm = numtype(sdm)==2 ? 0 : sdm

//Create _im waves in each datafolder for image plotting purposes
	duplicate/o $Nameofwave(diameter), $(Nameofwave(diameter)+"_im")
	wave diameter_im = $(Nameofwave(diameter)+"_im")
	duplicate/o $Nameofwave(dt), $(Nameofwave(dt)+"_im")
	wave dt_im = $(Nameofwave(dt)+"_im")
	insertpoints 0,1, diameter_im, dt_im
	//Use midpoints rather than edges
	diameter_im = diameter[p]-(diameter[p+1] - diameter[p])/2
	diameter_im[numpnts(diameter_im)-2] = diameter[numpnts(diameter)-1] - (diameter[numpnts(diameter)-1] - diameter[numpnts(diameter)-2])/2
	diameter_im[numpnts(diameter_im)-1] = diameter[numpnts(diameter)-1] + (diameter[numpnts(diameter)-1] - diameter[numpnts(diameter)-2])/2
	
//	dt_im = dt[p]-(dt[p+1] - dt[p])/2
//	dt_im[numpnts(dt_im)-2] = dt[numpnts(dt)-1] - (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
//	dt_im[numpnts(dt_im)-1] = dt[numpnts(dt)-1] + (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
	// do not use midpoints for time
	dt_im = dt[p]
	dt_im[numpnts(dt_im)-1] = dt[numpnts(dt)-1] + (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
	
//	diameter_im[0]=diameter_im[1]-(diameter_im[2]-diameter_im[1])
//	dt_im[0]=dt_im[1]-(dt_im[2]-dt_im[1])
	SetScale d 0,0,"dat", dt, dt_im
		
	//	dowindow/k Size_dist_image_plot
	//	execute("Size_dist_image_plot()")
	
	//	make/o/d/n=(dimsize(sdm,1)) $("recent_"+stringfromlist(0,S_fileName,"."))
	//	wave recent = $("recent_"+stringfromlist(0,S_fileName,"."))
	//	recent =sdm[dimsize(sdm,0)][p]
	//	dowindow/k Most_Recent_SMPS_dist_plot
	//	execute("Most_Recent_SMPS_dist("+nameofwave(recent)+","+nameofwave(diameter)+")") 
		killwaves smps0	
	endfor	
//	setdatafolder sdf
end

//************************************************************************************
function load_aps_files()
	// load a single APS file
	string extension = "txt"
	Variable refNum
	String message = "Select one or more files"
	String file_list
	String fileFilters = "Data Files (*."+extension+"):."+extension+";"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	file_list = S_fileName
	if (strlen(file_list) == 0)
		print "==================================="
		print "No APS Data Files Chosen for Loading.  Aborting."
		print "==================================="
		Close/A
		setdatafolder root:
		abort
	endif
	Close/A
	file_list = replacestring("\r",file_list,";")
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:aps	
	load_aps_row(file_list)
	setdatafolder root:aps
end

//************************************************************************************
function load_aps_directory()
	// load a directory of APS files
	newpath/o/q path_to_data	//, local_path_to_data
	string file_list = IndexedFile(path_to_data, -1, ".txt")
	if (cmpstr(file_list, "")==0)
		print "==================================="
		print "No APS *.txt data found in directory.  Aborting."
		print "==================================="
		setdatafolder root:
		abort
	endif
	//Enter *.txt filenames to remove from list (i.e. SKIP PROCESSING) using format in line below
//	file_list = RemoveFromList("LongDMA_calnex_SanDiego_11May2010_spl_400.txt", file_list)
	file_list = replacestring("\r",file_list,";")
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:aps	
//	return file_list	
	load_aps_row(file_list)
	setdatafolder root:aps
end

//************************************************************************************
function load_aps_row(file_list)
	// this function will load data from an APS file from a TSI APS
	// all data are loaded as a giant string wave, and then parsed accordingly
	string file_list 	//= initialize_smps()
	string sdf = getdatafolder(1)
	string file_full_path, file_name, file_name_str
	variable nf = itemsinlist(file_list)
	variable i_nf, nc, i_header, i
	for (i_nf=0;i_nf<nf;i_nf+=1)
		file_full_path = stringfromlist(i_nf,file_list)
		file_name = ParseFilePath(0, file_full_path, ":", 1, 0)
		file_name_str = ParseFilePath(3, file_full_path, ":", 0, 0)

		 if (strlen(file_name)>31)
//		 	file_name = file_name[0,30]
		 endif
		 
	//Test if the filename starts with number or not
		 if (numtype( str2num(file_name[0]) ) ==0)
		 	file_name_str = ("root:aps:x"+file_name_str)
		 else
		 	file_name_str = ("root:aps:"+file_name_str)
		 endif
		newdatafolder/o/s $file_name_str
		LoadWave/J/K=2/V={""," $",0,0}/L={0,0,0,0,1}/n=aps file_full_path
		wave/t aps0
		variable nr=numpnts(aps0)
	//Header INFO: Test if file contains header info, if so count header lines for appropriate wave, and Load Header Parameters
		if (numtype(str2num(aps0[0]))!=0)
			i_header=0
			for (i=0;i<nr;i+=1)
				if (numtype(str2num(aps0[i]))!=0)
					i_header+=1
				else
					break
				endif
			endfor
			//hdr = header info
			make/o/n=(i_header-1,2)/t  $("hdr_"+stringfromlist(0,S_fileName,"."))
			wave/t aps_header_info = $("hdr_"+stringfromlist(0,S_fileName,"."))
			aps_header_info[][0] = stringfromlist(0,aps0[p],"\t")
			aps_header_info[][1] = stringfromlist(1,aps0[p],"\t")
		endif			
		
	//Load Size Distribution Waves
		string d,t
		variable np=0,d1,d2,d3,t1,t2,t3
		//sn = sample number
		make/o/d/n=0 $("sn_"+stringfromlist(0,S_fileName,"."))
		wave sample_no = $("sn_"+stringfromlist(0,S_fileName,"."))
		make/o/d/n=0 $("dt_"+stringfromlist(0,S_fileName,"."))
		wave dt = $("dt_"+stringfromlist(0,S_fileName,"."))
		make/o/d/n=0 $("tdiff_"+stringfromlist(0,S_fileName,"."))
		wave tdiff = $("tdiff_"+stringfromlist(0,S_fileName,"."))
		for (i=0;i<nr;i+=1)												// Get Sample Number and DateTime
			if (numtype(str2num(aps0[i]))==0) 							// line starts with a number
				insertpoints np,1,dt,sample_no,tdiff
				sample_no[np]=str2num(stringfromlist(0,aps0[i],"\t"))
				d=stringfromlist(1,aps0[i],"\t")
				d1=str2num(stringfromlist(0,d,"/"))
				d2=str2num(stringfromlist(1,d,"/"))
				d3=str2num(stringfromlist(2,d,"/"))
				if(d3<2000)
					d3 += 2000 // in case year loads as YY instead of as YYYY
				endif
				t=stringfromlist(2,aps0[i],"\t")
				t1=str2num(stringfromlist(0,t,":"))
				t2=str2num(stringfromlist(1,t,":"))
				t3=str2num(stringfromlist(2,t,":"))
				dt[np]=date2secs(d3,d1,d2)+t1*3600+t2*60+t3
				tdiff[np] = dt[np] - dt[(np-1)]
				if (np==1)
					tdiff[0] = tdiff[np]
				endif
				np+=1
			endif
		endfor
		make/o/d/n=0 $("dia_"+stringfromlist(0,S_fileName,"."))
		wave diameter = $("dia_"+stringfromlist(0,S_fileName,"."))
		variable nd,j,k=0,l,m,namesline, start_sdm=0
		for (i=0;i<nr;i+=1) 											// loop over all lines
			if (stringmatch(stringfromlist(0,aps0[i],"\t"),"Sample #")) 		// look for line that starts with "Sample #"
				nd=itemsinlist(aps0[i],"\t")
				namesline = i
				for (j=0;j<nd;j+=1)										// loop over number of items in that line
					if (numtype(str2num(stringfromlist(j,aps0[i],"\t")))==0) 	// if the value is a number, assume it is a diameter; for APS ignore value with diameter as <0.XXX
						if(start_sdm==0)
							start_sdm=j								// index of first column with data
						endif
						insertpoints k,1,diameter						// build your wave
						diameter[k]=str2num(stringfromlist(j,aps0[i],"\t"))	// populate your wave
						k+=1
					endif
				endfor
				break
			endif
		endfor
		// get sample type (not necessary for SMPS, only APS)
		for(i=0;i<nd;i+=1)
			if(stringmatch(stringfromlist(i,aps0[namesline],"\t"),"Aerodynamic Diameter (mid)"))
				string sampletype = stringfromlist(i,aps0[namesline+1],"\t")
				insertpoints dimsize(aps_header_info,0),1,aps_header_info
				aps_header_info[dimsize(aps_header_info,0)][0] = "Units"
				aps_header_info[dimsize(aps_header_info,0)][1] = sampletype
			endif
		endfor
		//sdm = size distribution matrix
		make/o/d/n=(0,k) $("sdm_"+stringfromlist(0,S_fileName,"."))
		wave sdm = $("sdm_"+stringfromlist(0,S_fileName,"."))
		m=0
		l=0
		for (i=0;i<nr;i+=1)	
			if (numtype(str2num(aps0[i]))==0)
				nd=itemsinlist(aps0[i],"\t")
				insertpoints l,1,sdm
				for (j=start_sdm;j<k+start_sdm;j+=1)
					sdm[l][m]=str2num(stringfromlist(j,aps0[i],"\t"))
					m+=1
				endfor
				l+=1
				m=0
			endif
		endfor
		// turn NaN values into 0 values
		sdm = numtype(sdm)==2 ? 0 : sdm

//Create _im waves in each datafolder for image plotting purposes
	duplicate/o $Nameofwave(diameter), $(Nameofwave(diameter)+"_im")
	wave diameter_im = $(Nameofwave(diameter)+"_im")
	duplicate/o $Nameofwave(dt), $(Nameofwave(dt)+"_im")
	wave dt_im = $(Nameofwave(dt)+"_im")
	insertpoints 0,1, diameter_im, dt_im
	//Use midpoints rather than edges
	diameter_im = diameter[p]-(diameter[p+1] - diameter[p])/2
	diameter_im[numpnts(diameter_im)-2] = diameter[numpnts(diameter)-1] - (diameter[numpnts(diameter)-1] - diameter[numpnts(diameter)-2])/2
	diameter_im[numpnts(diameter_im)-1] = diameter[numpnts(diameter)-1] + (diameter[numpnts(diameter)-1] - diameter[numpnts(diameter)-2])/2
	
//	dt_im = dt[p]-(dt[p+1] - dt[p])/2
//	dt_im[numpnts(dt_im)-2] = dt[numpnts(dt)-1] - (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
//	dt_im[numpnts(dt_im)-1] = dt[numpnts(dt)-1] + (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
	// do not use midpoints for time
	dt_im = dt[p]
	dt_im[numpnts(dt_im)-1] = dt[numpnts(dt)-1] + (dt[numpnts(dt)-1] - dt[numpnts(dt)-2])/2
		
//	diameter_im[0]=diameter_im[1]-(diameter_im[2]-diameter_im[1])
//	dt_im[0]=dt_im[1]-(dt_im[2]-dt_im[1])
	SetScale d 0,0,"dat", dt, dt_im
		
	//	dowindow/k Size_dist_image_plot
	//	execute("Size_dist_image_plot()")
	
	//	make/o/d/n=(dimsize(sdm,1)) $("recent_"+stringfromlist(0,S_fileName,"."))
	//	wave recent = $("recent_"+stringfromlist(0,S_fileName,"."))
	//	recent =sdm[dimsize(sdm,0)][p]
	//	dowindow/k Most_Recent_aps_dist_plot
	//	execute("Most_Recent_aps_dist("+nameofwave(recent)+","+nameofwave(diameter)+")") 
		killwaves aps0	
	endfor	
//	setdatafolder sdf
end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions to load time-series of CPC number concentrations determined
// using a TSI CPC and exported from the AIM program.
// Text files must be exported as *.txt, columns, time stamp
// Export only the concentration, and comma delimited
function load_cpc_files()
// load a selection of CPC files
	string extension = "txt"
	Variable refNum
	String message = "Select one or more files"
	String file_list
	String fileFilters = "Data Files (*."+extension+"):."+extension+";"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	file_list = S_fileName
	if (strlen(file_list) == 0)
		print "==================================="
		print "No Data Files Chosen for Loading.  Aborting."
		print "==================================="
		Close/A
		setdatafolder root:
		abort
	endif
	Close/A
	file_list = replacestring("\r",file_list,";")
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:cpc	
//	load_cpc_row(file_list)
	load_cpc_row_comma(file_list)
	setdatafolder root:cpc
end

//************************************************************************************
function load_cpc_directory()
	newpath/o/q path_to_data	//, local_path_to_data
	string file_list = IndexedFile(path_to_data, -1, ".txt")
	if (cmpstr(file_list, "")==0)
		print "==================================="
		print "No cpc *.txt data found in directory.  Aborting."
		print "==================================="
		setdatafolder root:
		abort
	endif
	//Enter *.txt filenames to remove from list (i.e. SKIP PROCESSING) using format in line below
	file_list = sortlist(file_list,";",16)
	newdatafolder/o/s root:cpc	
	load_cpc_row(file_list)
	setdatafolder root:cpc
end

//************************************************************************************
function load_cpc_row(file_list)
// Parse the selection of CPC files into time series. 
// NaN's are deleted (sometimes these hang around due to problems with loading)
	string file_list 	
	string sdf = getdatafolder(1)
	string file_full_path, file_name, file_name_str
	variable nf = itemsinlist(file_list)
	variable i_nf, nc, i_header, i, j
	
	variable dateline = 4	// line number on which date information for each scan lives
	
	for (i_nf=0;i_nf<nf;i_nf+=1)	// loop over all files loaded
		file_full_path = stringfromlist(i_nf,file_list)
		file_name = ParseFilePath(0, file_full_path, ":", 1, 0)
		file_name_str = ParseFilePath(3, file_full_path, ":", 0, 0)

		 if (strlen(file_name)>31)
//		 	file_name = file_name[0,30]
		 endif
		 
	//Test if the filename starts with number or not
		 if (numtype( str2num(file_name[0]) ) ==0)
		 	file_name_str = ("root:cpc:x"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 else
		 	file_name_str = ("root:cpc:"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 endif
		newdatafolder/o/s $file_name_str
		LoadWave/J/K=2/V={""," $",0,0}/L={0,0,0,0,1}/n=cpc file_full_path
//		LoadWave/J/K=2/V={","," $",0,0}/L={0,0,0,0,1}/n=cpc file_full_path

		wave/t cpc0
		variable nr=numpnts(cpc0)
	//Header INFO: Test if file contains header info, if so count header lines for appropriate wave, and Load Header Parameters
		if (numtype(str2num(cpc0[0]))!=0)
			i_header=0
			for (i=0;i<nr;i+=1)
				if (numtype(str2num(cpc0[i]))!=0)
					i_header+=1
				else
					break
				endif
			endfor
			//hdr = header info
			make/o/n=(i_header-1,2)/t  $("hdr_"+stringfromlist(0,S_fileName,"."))
			wave/t cpc_header_info = $("hdr_"+stringfromlist(0,S_fileName,"."))
			cpc_header_info[][0] = stringfromlist(0,cpc0[p],"\t")
			cpc_header_info[][1] = stringfromlist(1,cpc0[p],"\t")
		endif			

	variable n_samples = itemsinlist(cpc0[i_header],"\t")/2 // number of scans
	string currentstr, firstcharstr
	variable i_stop // index for last point
	for(i=i_header;i<nr;i+=1) // get index for last data point
		currentstr = cpc0[i]
		firstcharstr = currentstr[0]
		if(numtype(str2num(firstcharstr))==2)
			i_stop = i - 1
			break
		endif
	endfor
	
	// Make a bunch of silly waves to deal with text parsing. Overly complicated, ultimately, but it works
	make/o/d/n=((i_stop-i_header)*n_samples) $("hr_"+stringfromlist(0,S_fileName,"."))
	wave hr_single = $("hr_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("min_"+stringfromlist(0,S_fileName,"."))
	wave min_single = $("min_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("sec_"+stringfromlist(0,S_fileName,"."))
	wave sec_single = $("sec_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("mo_"+stringfromlist(0,S_fileName,"."))
	wave mo_single = $("mo_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("day_"+stringfromlist(0,S_fileName,"."))
	wave day_single = $("day_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("yr_"+stringfromlist(0,S_fileName,"."))
	wave yr_single = $("yr_"+stringfromlist(0,S_fileName,"."))
	make/o/t/n=((i_stop-i_header)*n_samples) $("tstr_"+stringfromlist(0,S_fileName,".")) = ""
	wave/t tstr_single = $("tstr_"+stringfromlist(0,S_fileName,"."))
	// Make your final time and particle number waves
	make/o/d/n=((i_stop-i_header)*n_samples) $("ts_"+stringfromlist(0,S_fileName,"."))
	wave t_single = $("ts_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("Nps_"+stringfromlist(0,S_fileName,".")) = nan
	wave Np_single = $("Nps_"+stringfromlist(0,S_fileName,"."))
	
	variable hr, minute, second, yr
	string currentdate
	variable counter
//	// Loop to populate waves for individual scans
//	for(i=0;i<n_samples;i+=1) // num samples
//		currentdate = stringfromlist(i*2+1,cpc0[dateline],"\t")
//		for(j=0;j<(i_stop-i_header);j+=1) // num rows
//			counter = j+i*(i_stop-i_header)
//			currentstr = (stringfromlist(i*2,cpc0[j+i_header],"\t"))
//			tstr_single[j+i*(i_stop-i_header)] = currentstr
//			hr_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(0,currentstr,":"))
//			min_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(1,currentstr,":"))
//			sec_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(2,currentstr,":"))
//			mo_single[counter] = str2num(stringfromlist(0,currentdate,"/"))
//			day_single[counter] = str2num(stringfromlist(1,currentdate,"/"))
//			yr = str2num(stringfromlist(2,currentdate,"/"))
//			yr_single[counter] = yr < 2000 ? yr+2000 : yr
//			t_single[counter] = date2secs(yr_single[counter],mo_single[counter],day_single[counter])+60*60*hr_single[counter]+60*min_single[counter]+sec_single[counter]
//			Np_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(i*2+1,cpc0[j+i_header],"\t"))
//		endfor
//	endfor

	variable t1, t2
	variable dt = 0
	variable i_newday
	variable da, mo
	
	// Loop to populate waves for individual scans
	for(i=0;i<n_samples;i+=1) // num samples
		currentdate = stringfromlist(i*2+1,cpc0[dateline],"\t")
		dt = 0
		for(j=0;j<(i_stop-i_header);j+=1) // num rows
			counter = j+i*(i_stop-i_header)
			currentstr = (stringfromlist(i*2,cpc0[j+i_header],"\t"))
			tstr_single[j+i*(i_stop-i_header)] = currentstr
			hr_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(0,currentstr,":"))
			min_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(1,currentstr,":"))
			sec_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(2,currentstr,":"))
			
			if(j==0)
				t1 = 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
				mo = str2num(stringfromlist(0,currentdate,"/"))
				da = str2num(stringfromlist(1,currentdate,"/"))
				yr = str2num(stringfromlist(2,currentdate,"/"))
				yr = yr < 2000 ? yr + 2000 : yr

			elseif(j==1)
				t2 = 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
				dt = t2-t1
			endif
			t_single[counter] = date2secs(yr,mo,da) + t1 + dt*j // 60*60*hr_single[counter]+60*min_single[counter]+sec_single[counter]
			Np_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(i*2+1,cpc0[j+i_header],"\t"))
		endfor
	endfor
	
	setscale/p d, 0, 0, "dat", t_single
	
	killwaves tstr_single, hr_single, min_single, sec_single, mo_single, day_single, yr_single
	endfor	
////	setdatafolder sdf
end


//************************************************************************************
function load_cpc_row_comma(file_list)
// Parse the selection of CPC files into time series. 
// NaN's are deleted (sometimes these hang around due to problems with loading)
	string file_list 	
	string sdf = getdatafolder(1)
	string file_full_path, file_name, file_name_str
	variable nf = itemsinlist(file_list)
	variable i_nf, nc, i_header, i, j
	
	variable dateline = 4	// line number on which date information for each scan lives
	variable startline = 5
	variable lengthline = 6
	variable intervalline = 7
	
	for (i_nf=0;i_nf<nf;i_nf+=1)	// loop over all files loaded
		file_full_path = stringfromlist(i_nf,file_list)
		file_name = ParseFilePath(0, file_full_path, ":", 1, 0)
		file_name_str = ParseFilePath(3, file_full_path, ":", 0, 0)

		 if (strlen(file_name)>31)
//		 	file_name = file_name[0,30]
		 endif
		 
	//Test if the filename starts with number or not
		 if (numtype( str2num(file_name[0]) ) ==0)
		 	file_name_str = ("root:cpc:x"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 else
		 	file_name_str = ("root:cpc:"+file_name_str)	//file_name[0,strlen(file_name)-5])
		 endif
		newdatafolder/o/s $file_name_str
		LoadWave/J/K=2/V={""," $",0,0}/L={0,0,0,0,1}/n=cpc file_full_path
//		LoadWave/J/K=2/V={","," $",0,0}/L={0,0,0,0,1}/n=cpc file_full_path

		wave/t cpc0
		variable nr=numpnts(cpc0)
	//Header INFO: Test if file contains header info, if so count header lines for appropriate wave, and Load Header Parameters
		if (numtype(str2num(cpc0[0]))!=0)
			i_header=0
			for (i=0;i<nr;i+=1)
				if (numtype(str2num(cpc0[i]))!=0)
					i_header+=1
				else
					break
				endif
			endfor
			//hdr = header info
			make/o/n=(i_header-1,2)/t  $("hdr_"+stringfromlist(0,S_fileName,"."))
			wave/t cpc_header_info = $("hdr_"+stringfromlist(0,S_fileName,"."))
			cpc_header_info[][0] = stringfromlist(0,cpc0[p],",")
			cpc_header_info[][1] = stringfromlist(1,cpc0[p],",")
		endif			

	variable n_samples = itemsinlist(cpc0[i_header],",")/2 // number of scans
	string currentstr, firstcharstr
	variable i_stop // index for last point
	for(i=i_header;i<nr;i+=1) // get index for last data point
		currentstr = cpc0[i]
		firstcharstr = currentstr[0]
		if(stringmatch(firstcharstr,"C"))
//		if(numtype(str2num(firstcharstr))==2)
			i_stop = i-2 // i - 1
			break
		endif
	endfor
	
	// Make a bunch of silly waves to deal with text parsing. Overly complicated, ultimately, but it works
	make/o/d/n=((i_stop-i_header)*n_samples) $("hr_"+stringfromlist(0,S_fileName,"."))
	wave hr_single = $("hr_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("min_"+stringfromlist(0,S_fileName,"."))
	wave min_single = $("min_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("sec_"+stringfromlist(0,S_fileName,"."))
	wave sec_single = $("sec_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("mo_"+stringfromlist(0,S_fileName,"."))
	wave mo_single = $("mo_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("day_"+stringfromlist(0,S_fileName,"."))
	wave day_single = $("day_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=((i_stop-i_header)*n_samples) $("yr_"+stringfromlist(0,S_fileName,"."))
	wave yr_single = $("yr_"+stringfromlist(0,S_fileName,"."))
	make/o/t/n=((i_stop-i_header)*n_samples) $("tstr_"+stringfromlist(0,S_fileName,".")) = ""
	wave/t tstr_single = $("tstr_"+stringfromlist(0,S_fileName,"."))
	// Make your final time and particle number waves
//	make/o/d/n=((i_stop-i_header)*n_samples) $("ts_"+stringfromlist(0,S_fileName,"."))
//	wave t_single = $("ts_"+stringfromlist(0,S_fileName,"."))
//	make/o/d/n=((i_stop-i_header)*n_samples) $("Nps_"+stringfromlist(0,S_fileName,".")) = nan
//	wave Np_single = $("Nps_"+stringfromlist(0,S_fileName,"."))
	// Make your final time and particle number waves
	make/o/d/n=(0) $("ts_"+stringfromlist(0,S_fileName,"."))
	wave t_single = $("ts_"+stringfromlist(0,S_fileName,"."))
	make/o/d/n=(0) $("Nps_"+stringfromlist(0,S_fileName,".")) = nan
	wave Np_single = $("Nps_"+stringfromlist(0,S_fileName,"."))
	
	variable hr, minute, second, yr
	string currentdate
	variable counter
	string currentstr_new, first
	variable stringlength
	
	// remove unnecessary commas that muck things up
	for(j=0;j<(i_stop-i_header);j+=1) // num rows
		counter = j+i*(i_stop-i_header)
		currentstr = cpc0[j+i_header]
		currentstr_new = ""
		first = currentstr[0,3]
		if(stringmatch(first,",,,,"))
			currentstr_new = currentstr[4,strlen(currentstr)-4]
			currentstr = currentstr_new
		endif
		cpc0[j+i_header] = currentstr
	endfor

	variable t1, t2
	variable dt = 0
	variable i_newday
	variable da, mo
	
//	// Loop to populate waves for individual scans
//	for(i=0;i<n_samples;i+=1) // num samples
//		currentdate = stringfromlist(i*2+1,cpc0[dateline],",")
//		dt = 0
//		for(j=0;j<(i_stop-i_header);j+=1) // num rows
//			counter = j+i*(i_stop-i_header)
//			currentstr = (stringfromlist(i*2,cpc0[j+i_header],","))
//			tstr_single[j+i*(i_stop-i_header)] = currentstr
//			hr_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(0,currentstr,":"))
//			min_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(1,currentstr,":"))
//			sec_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(2,currentstr,":"))
//			
//			if(j==0)
//				t1 = 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
//				mo = str2num(stringfromlist(0,currentdate,"/"))
//				da = str2num(stringfromlist(1,currentdate,"/"))
//				yr = str2num(stringfromlist(2,currentdate,"/"))
//				yr = yr < 2000 ? yr + 2000 : yr
//
//			elseif(j==1)
//				t2 = 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
//				dt = t2-t1
//			endif
//			t_single[counter] = date2secs(yr,mo,da) + t1 + dt*j // 60*60*hr_single[counter]+60*min_single[counter]+sec_single[counter]
//			Np_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(i*2+1,cpc0[j+i_header],","))
//		endfor
//	endfor
	
	string currentstarttime
	string currentlength
	variable currentinterval
	variable hr_length, min_length, sec_length
	variable hr_start, min_start, sec_start
	variable i_length
	variable current_length
	
		// Loop to populate waves for individual scans
	for(i=0;i<n_samples;i+=1) // num samples
		currentdate = stringfromlist(i*2+1,cpc0[dateline],",")
		currentstarttime = stringfromlist(i*2+1,cpc0[dateline+1],",")
		currentlength = stringfromlist(i*2+1,cpc0[dateline+2],",")
//		print currentlength
		dt = str2num(stringfromlist(i*2+1,cpc0[dateline+3],","))
		hr_start = str2num(stringfromlist(0,currentstarttime,":"))
		min_start = str2num(stringfromlist(1,currentstarttime,":"))
		sec_start = str2num(stringfromlist(2,currentstarttime,":"))
		hr_length = str2num(stringfromlist(0,currentlength,":"))
		min_length = str2num(stringfromlist(1,currentlength,":"))
		sec_length = str2num(stringfromlist(2,currentlength,":"))
		if(numtype(sec_length)==2) // is nan
			sec_length = min_length
			min_length = hr_length
			hr_length = 0
		elseif(numtype(min_length)==2) // is nan
			sec_length = hr_length
			min_length = 0
			hr_length = 0
		endif
		i_length = (hr_length*60*60+min_length*60+sec_length)/dt
//		print i_length
		
		current_length = numpnts(t_single)
		redimension/n=(current_length+i_length) t_single, Np_single
		
//		dt = 0
		for(j=0;j<i_length;j+=1) // num rows
//			counter = j+i*(i_stop-i_header)
			currentstr = (stringfromlist(i*2,cpc0[j+i_header],","))
//			tstr_single[j+i*(i_stop-i_header)] = currentstr
//			hr_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(0,currentstr,":"))
//			min_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(1,currentstr,":"))
//			sec_single[j+i*(i_stop-i_header)] = str2num(stringfromlist(2,currentstr,":"))
			t1 = 3600*hr_start + 60*min_start + sec_start
			if(j==0)
				 // 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
				mo = str2num(stringfromlist(0,currentdate,"/"))
				da = str2num(stringfromlist(1,currentdate,"/"))
				yr = str2num(stringfromlist(2,currentdate,"/"))
				yr = yr < 2000 ? yr + 2000 : yr

			elseif(j==1)
//				t2 = 3600*hr_single[j+i*(i_stop-i_header)] + 60*min_single[j+i*(i_stop-i_header)] + sec_single[j+i*(i_stop-i_header)]
//				dt = t2-t1
			endif
			t_single[current_length+j] = date2secs(yr,mo,da) + t1 + dt*j // 60*60*hr_single[counter]+60*min_single[counter]+sec_single[counter]
			Np_single[current_length+j] = str2num(stringfromlist(i*2+1,cpc0[j+i_header],","))
		endfor
	endfor
	
	setscale/p d, 0, 0, "dat", t_single
	
	killwaves tstr_single, hr_single, min_single, sec_single, mo_single, day_single, yr_single
	endfor	
////	setdatafolder sdf
end

//************************************************************************************
Function MergeCPC_tseries()

	string df_CPC = "root:CPC:"
	setdatafolder $df_CPC
	string df_CPC_Merged = "root:CPC:Merged"
	if(datafolderexists(df_CPC_Merged)==0)
		newdatafolder $df_CPC_Merged
	endif
	
	variable ndf = countobjects(":",4)	//directories with files
	string datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("Merged") datafolder name from list
	if (stringmatch(datafolder_list,"*Merged*"))
		datafolder_list = RemoveFromList("Merged",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	
	
	setdatafolder df_CPC_merged

	make/o/n=(ndf) wv_npnts = nan
	
	variable V_npntstot 
	variable V_dia
	
	string cdf_str
	string wvstr_Np = "Nps_"
	string wvstr_time = "ts_"
	string current_str
	
	variable i, j
		
	string dt_list = "", Np_list=""
	
	for(i=0;i<ndf;i+=1)
		cdf_str = df_CPC + stringfromlist(i,datafolder_list,",")//[i]
		setdatafolder $cdf_str
		
		dt_list += df_CPC+stringfromlist(i,datafolder_list,",")+":"
		dt_list += wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		
		Np_list += df_CPC+stringfromlist(i,datafolder_list,",")+":"
		Np_list += wvstr_Np+stringfromlist(i,datafolder_list,",")+";"

	endfor
	
	setdatafolder df_CPC_Merged	
	concatenate/o dt_list, t_CPC_temp
	concatenate/o Np_list, Np_CPC_temp
	
	extract/o t_cpc_temp, t_CPC, numtype(t_cpc_temp)!=2
	extract/o Np_cpc_temp, Np_CPC, numtype(t_cpc_temp)!=2
	
	setscale/p d 0,1, "dat", t_CPC
		
	killwaves/z wv_npnts, t_cpc_temp, Np_cpc_temp
	
	setdatafolder df_CPC_Merged	
End

///// END OF CPC LOADING AND MERGING FILES /////

//************************************************************************************
function process_smps_data()
	// Tim Onasch's processing stuff
//&	setdatafolder root:smps
	string cdf = getdatafolder(1)	//Get Current Data Folder!

//TSI SMPS accounts for CPC flow dilution by comparing the input flowrate with the nominal instrument flow rate
	
//Concatenate all files together
	//find and concatenate all waves
	variable i,j
	variable ndf = countobjects(":",4)	//directories with files
	string datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("avg_dist") datafolder name from list
	if (stringmatch(datafolder_list,"*avg_dist*"))
		datafolder_list = RemoveFromList("avg_dist",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	
		
	//Lists of waves to be concatenated
	string cat_wave_list = "date_time_list;sample_no_list;time_diff_list;"
	string cat_wavename_list = "smps_date_time;smps_sample_no;smps_sample_time;"
	string/g date_time_list = ""
	string/g sample_no_list = ""
	string/g time_diff_list = ""
	//Lists of matrices to be concatenated
	string matrix_list = ""	
	
	for (i=0;i<ndf;i+=1)
		datafolder_str=stringfromlist(i,datafolder_list,",")	//note the string is separated by commas...	//GetIndexedObjName(":",4,j)
//		datafolder_str = GetIndexedObjName(":",4,i)

	//Check to see if "x" has been added to datafolder due to start of file name being a number...
		if (stringmatch(datafolder_str,"x*")==1)
			file_name_str = datafolder_str[1,strlen(datafolder_str)]
		else
			file_name_str = datafolder_str
		endif
//		if (stringmatch(datafolder_str,"avg_dist"))
//			ndf -= 1
//			continue
//		endif
		
	//Assume that all data uses same diameter wave, if NOT, then complain and quit...
//&		wave/z smps_diameter = $("root:smps:smps_diameter")
		wave/z smps_diameter = $(cdf+"smps_diameter")
		if (!waveexists(smps_diameter))		//(i==0)
//&			duplicate/o $("root:smps:"+datafolder_str+":dia_"+file_name_str), $("root:smps:smps_diameter")
//&			duplicate/o $("root:smps:"+datafolder_str+":dia_"+file_name_str+"_im"), $("root:smps:smps_diameter_im")
			duplicate/o $(cdf+datafolder_str+":dia_"+file_name_str), $("smps_diameter")
			duplicate/o $(cdf+datafolder_str+":dia_"+file_name_str+"_im"), $("smps_diameter_im")
		else
			//check if the diameter waves have changed between the data!!!!
//&			wave smps_diameter = $("root:smps:smps_diameter")
			wave smps_diameter = $(cdf+"smps_diameter")
//&			wave source_dia = $("root:smps:"+datafolder_str+":dia_"+file_name_str)
			wave source_dia = $(cdf+datafolder_str+":dia_"+file_name_str)
			if (!equalWaves(source_dia, smps_diameter, 1))
				print "Diameter waves are NOT the same!  Must process separately.  Quiting."
				abort
			endif
		endif

	//Generate lists of waves to be concatenated
//&		date_time_list += "root:smps:"+datafolder_str+":dt_"+file_name_str+";"
//&		sample_no_list += "root:smps:"+datafolder_str+":sn_"+file_name_str+";"
//&		time_diff_list += "root:smps:"+datafolder_str+":tdiff_"+file_name_str+";"
		date_time_list += cdf+datafolder_str+":dt_"+file_name_str+";"
		sample_no_list += cdf+datafolder_str+":sn_"+file_name_str+";"
		time_diff_list += cdf+datafolder_str+":tdiff_"+file_name_str+";"

	//Generate lists of matrices to be concatenated
//&		matrix_list += "root:smps:"+datafolder_str+":sdm_"+file_name_str+";"		
		matrix_list += cdf+datafolder_str+":sdm_"+file_name_str+";"		
	endfor

//Generate concatenated waves
	variable break_time = 2	//Time inserted between two different SMPS files to provide NaN breaks in visual data (seconds)
	variable num_wl = itemsinlist(cat_wave_list)
	string str_to_concatenate
	string/g source_list
	for (i=0;i<num_wl;i+=1)
		SVAR source_list = $stringfromlist(i,cat_wave_list)		//source_list = string containing a string reference...
	 	str_to_concatenate = ""
		for (j=0;j<ndf;j+=1)
			wave smps1 = $stringfromlist(j,source_list)
			wave/z smps2 = $stringfromlist((j+1),source_list)
			duplicate/o smps1, $(nameofwave(smps1)+"_add")
			wave smps3 = $(nameofwave(smps1)+"_add")
			if (j<(itemsinlist(source_list)-1))
				insertpoints numpnts(smps1),2,smps3
				smps3[numpnts(smps1)] = smps1[numpnts(smps1)-1]+break_time	
				smps3[numpnts(smps1)+1] = smps2[0]-break_time
			endif
			str_to_concatenate += nameofwave(smps3)+";"
		endfor
//&		concatenate /NP/O str_to_concatenate,$("root:smps:"+stringfromlist(i,cat_wavename_list))
		concatenate /NP/O str_to_concatenate,$(cdf+stringfromlist(i,cat_wavename_list))
		kill_list_of_waves_smps(wavelist("*_add",";",""))
	 endfor

//Make SMPS_start_time
//&	wave smps_date_time = $("root:smps:smps_date_time")
//&	wave smps_sample_time = $("root:smps:smps_sample_time")
//&	duplicate/o smps_date_time, $("root:smps:smps_start_time")
//&	wave smps_start_time = $("root:smps:smps_start_time")
	wave smps_date_time = $(cdf+"smps_date_time")
	wave smps_sample_time = $(cdf+"smps_sample_time")
	duplicate/o smps_date_time, $(cdf+"smps_start_time")
	wave smps_start_time = $(cdf+"smps_start_time")
	
	variable diff = median_smps(smps_sample_time)
	for (i=0;i<numpnts(sample_time);i+=1)
		if ((smps_sample_time[i] > (smps_sample_time[i+1]+diff)) && (smps_sample_time[i] > (smps_sample_time[i-1]+diff)))
			smps_sample_time[i] = smps_sample_time[i+1]
		endif
	endfor
	smps_start_time = smps_date_time - smps_sample_time
	
//Generate concatenated N matrix
	matrix_list = sortlist(matrix_list,";",16)
	variable num1=0,num2=0
//&	wave smps_diameter = $("root:smps:smps_diameter")
//&	make/o/n=(0,numpnts(smps_diameter)) $("root:smps:smps_Nmatrix")	//+num2str(i))
//&	wave smps_Nmatrix = $("root:smps:smps_Nmatrix")	//+num2str(i))		
	wave smps_diameter = $(cdf+"smps_diameter")
	make/o/n=(0,numpnts(smps_diameter)) $(cdf+"smps_Nmatrix")	//+num2str(i))
	wave smps_Nmatrix = $(cdf+"smps_Nmatrix")	//+num2str(i))		

	for (i=0;i<ndf;i+=1)
		wave smps1 = $stringfromlist(i,matrix_list)
		num2 += dimsize(smps1,0)
		insertpoints/m=0 num1,(num2-num1),smps_Nmatrix
		smps_Nmatrix[num1,num2-1][] = smps1[p-num1][q]
		if (i<(itemsinlist(matrix_list)-1))
			insertpoints/m=0 num2,2,smps_Nmatrix
			smps_Nmatrix[num2][] = NaN
			smps_Nmatrix[num2+1][] = NaN
			num2+=2
		endif
		num1 = num2		
	endfor

//Generate S matrix - surface area matrix
//&	duplicate/o smps_Nmatrix,$("root:smps:smps_Smatrix")	//+num2str(i))
//&	wave smps_Smatrix = $("root:smps:smps_Smatrix")		//+num2str(i))
//&	wave smps_diameter = $("root:smps:smps_diameter")		//+num2str(i))
	duplicate/o smps_Nmatrix,$(cdf+"smps_Smatrix")	//+num2str(i))
	wave smps_Smatrix = $(cdf+"smps_Smatrix")		//+num2str(i))
	wave smps_diameter = $(cdf+"smps_diameter")		//+num2str(i))
	smps_Smatrix = pi * smps_diameter[q]^2 * smps_Nmatrix[p][q] / 10^14 		//units of cm2/cm3/dlogDp x 1/density (cm3/g) x dlogDp = cm2/g or mass specific surface area  //may wish to keep units of nm2/cm3/dlogDp		-

//Generate V matrix - volume matrix
//&	duplicate/o smps_Nmatrix,$("root:smps:smps_Vmatrix")	//+num2str(i))
//&	wave smps_Vmatrix = $("root:smps:smps_Vmatrix")		//+num2str(i))
//&	wave smps_diameter = $("root:smps:smps_diameter")		//+num2str(i))
	duplicate/o smps_Nmatrix,$(cdf+"smps_Vmatrix")	//+num2str(i))
	wave smps_Vmatrix = $(cdf+"smps_Vmatrix")		//+num2str(i))
	wave smps_diameter = $(cdf+"smps_diameter")		//+num2str(i))
	smps_Vmatrix = pi/6 * smps_diameter[q]^3 * smps_Nmatrix[p][q] / 10^9		//units of um3/cm3/dlogDp x density (1g/cm3) = ug/m3/dlogDp

//Calculate modes and means and make root:smps:*_im waves
	calculate_modes_and_means()

	KillStrings/A/Z

//Create SMPS matrix plot
//	dowindow/k SMPS_matrix_plot
//	execute("SMPS_matrix_plot()")
end

//************************************************************************************
function calculate_modes_and_means()
	// calculate modes and means from size distribution matrices
	variable i,j
//&	setdatafolder root:smps
	
//Integrate data and calculate mode and geometric mean diameters
	wave Nmatrix = smps_Nmatrix
	wave Smatrix = smps_Smatrix
	wave Vmatrix = smps_Vmatrix
	wave date_time = smps_date_time
	wave sample_time = smps_sample_time
	wave diameter = smps_diameter
	wave diameter_im = smps_diameter_im
	variable np=dimsize(smps_Nmatrix,0)
	variable nq=dimsize(smps_Nmatrix,1)
	//Create DlogDp wave
	make/o/n=(nq)/d $("smps_dlogDp")
	wave dlogDp = $("smps_dlogDp")
	dlogDp=log(diameter_im[p+1])-log(diameter_im[p])	//Now using new smps_diameter_im wave rather than smps_diameter
//	dlogDp[0]=dlogDp[1]	
	//Number mode info	
	make/o/d/n=(np) $("smps_integrated_N")
	wave integrated_N = $("smps_integrated_N")
	make/o/d/n=(np) $("smps_Nmode_dia")
	wave Nmode_dia = $("smps_Nmode_dia")
	make/o/d/n=(np) $("smps_NGeomean_dia")
	wave Ngeomean_dia = $("smps_NGeomean_dia")
	make/o/d/n=(np) $("smps_NGeoStdev")
	wave NGeoStdev = $("smps_NGeoStdev")
//	make/o/d/n=(np) $("smps_number_per_dia")
//	wave number_per_dia = $("smps_number_per_dia")
//	integrated_N=0
//	Ngeomean_dia = 0
//	number_per_dia = 0
	//Surface area mode info	
	make/o/d/n=(np) $("smps_integrated_S")
	wave integrated_S = $("smps_integrated_S")
	make/o/d/n=(np) $("smps_Smode_dia")
	wave Smode_dia = $("smps_Smode_dia")
	make/o/d/n=(np) $("smps_SGeomean_dia")
	wave Sgeomean_dia = $("smps_SGeomean_dia")
	make/o/d/n=(np) $("smps_SGeoStdev")
	wave SGeoStdev = $("smps_SGeoStdev")
//	make/o/d/n=(np) $("smps_suface_area_per_dia")
//	wave suface_area_per_dia = $("smps_suface_area_per_dia")
//	integrated_S=0
//	Sgeomean_dia = 0
//	suface_area_per_dia = 0
	//Volume mode info 	
	make/o/d/n=(np) $("smps_integrated_V")
	wave integrated_V = $("smps_integrated_V")
	make/o/n=(np) $("smps_Vmode_dia")
	wave Vmode_dia = $("smps_Vmode_dia")
	make/o/d/n=(np) $("smps_VGeomean_dia")
	wave Vgeomean_dia = $("smps_VGeomean_dia")
	make/o/d/n=(np) $("smps_VGeoStdev")
	wave VGeoStdev = $("smps_VGeoStdev")
//	make/o/d/n=(np) $("smps_volume_per_dia")
//	wave volume_per_dia = $("smps_volume_per_dia")
//	integrated_V = 0
//	Vgeomean_dia = 0
//	volume_per_dia = 0
	make/o/n=(nq) dist_N, dist_V, dist_S
	for (i=0;i<np;i+=1)
		dist_N = Nmatrix[i][p]
//		wavestats/q dist_N
//		Nmode_dia[i] = diameter[V_maxloc]
		wave/z N_dist_stats = calc_lognormal_statistics2(dist_N,diameter)
		integrated_N[i] = N_dist_stats[0]
		Nmode_dia[i] = N_dist_stats[1]
		Ngeomean_dia[i] = N_dist_stats[2]
		NGeoStdev[i] = N_dist_stats[3]

		dist_S = Smatrix[i][p]
//		wavestats/q dist_S
//		Smode_dia[i] = diameter[V_maxloc]
		wave/z S_dist_stats = calc_lognormal_statistics2(dist_S,diameter)
		integrated_S[i] = S_dist_stats[0]
		Smode_dia[i] = S_dist_stats[1]
		Sgeomean_dia[i] = S_dist_stats[2]
		SGeoStdev[i] = S_dist_stats[3]

		dist_V = Vmatrix[i][p]
//		wavestats/q dist_V
//		Vmode_dia[i] = diameter[V_maxloc]
		wave/z V_dist_stats = calc_lognormal_statistics2(dist_V,diameter)
		integrated_V[i] = V_dist_stats[0]
		Vmode_dia[i] = V_dist_stats[1]
		Vgeomean_dia[i] = V_dist_stats[2]
		VGeoStdev[i] = V_dist_stats[3]

//		for (j=0;j<nq;j+=1)
//			Ngeomean_dia[i] += dist_N[j]*dlogDp[j]*log(diameter[j])
//			integrated_N[i] += Nmatrix[i][j]*dlogDp[j]
//			number_per_dia[i] += dist_N[j]*dlogDp[j]
//
//			Sgeomean_dia[i] += dist_S[j]*dlogDp[j]*log(diameter[j])
//			integrated_S[i] += Smatrix[i][j]*dlogDp[j]
//			suface_area_per_dia[i] += dist_S[j]*dlogDp[j]
//
//			Vgeomean_dia[i] += dist_V[j]*dlogDp[j]*log(diameter[j])
//			integrated_V[i] += Vmatrix[i][j]*dlogDp[j]
//			volume_per_dia[i] += dist_V[j]*dlogDp[j]
//		endfor
//		Ngeomean_dia[i] = 10^(Ngeomean_dia[i]/number_per_dia[i])
//		Sgeomean_dia[i] = 10^(Sgeomean_dia[i]/suface_area_per_dia[i])
//		Vgeomean_dia[i] = 10^(Vgeomean_dia[i]/volume_per_dia[i])
	endfor
	killwaves dist_N,dist_S,dist_V	//, number_per_dia, volume_per_dia, suface_area_per_dia
	
//Create _im waves in datafolder root:smps: for image plotting purposes
//	duplicate/o $Nameofwave(diameter), $(Nameofwave(diameter)+"_im")
//	wave diameter_im = $(Nameofwave(diameter)+"_im")
	duplicate/o $Nameofwave(date_time), $(Nameofwave(date_time)+"_im")
	wave date_time_im = $(Nameofwave(date_time)+"_im")
	insertpoints 0,1, date_time_im		//diameter_im, 
//	diameter_im[0]=diameter_im[1]-(diameter_im[2]-diameter_im[1])

	//Use mid-points rather than off-set points and use the sample time to determine the time intervals, rather than just total time
	date_time_im = date_time[p] - sample_time[p]/2
//	date_time_im[numpnts(date_time_im)-2] = date_time[numpnts(date_time)-1] - sample_time[numpnts(sample_time)-1]/2
	date_time_im[numpnts(date_time_im)-1] = date_time[numpnts(date_time)-1] + sample_time[numpnts(sample_time)-1]/2
	
//	date_time_im[0]=date_time_im[1]-(date_time_im[2]-date_time_im[1])
//&	SetScale d 0,0,"dat", root:SMPS:smps_date_time,root:SMPS:smps_date_time_im
	SetScale d 0,0,"dat", smps_date_time,smps_date_time_im
end

//************************************************************************************
//SMPS_matrix_plot() is now just a macro model for the fxn SMPS_matrix_plot2()
//	Need to use the fxn form to deal with different datafolders!
// Tim Onasch's function
Window SMPS_matrix_plot() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(24.75,51.5,669,428) as "SMPS_matrix_plot"
	AppendImage root:analysis_old:smps_Nmatrix vs {smps_date_time_im,smps_diameter_im}
	ModifyImage smps_Nmatrix ctab= {*,50000,Rainbow,1}
	AppendImage/L=vol smps_Vmatrix vs {smps_date_time_im,smps_diameter_im}
	ModifyImage smps_Vmatrix ctab= {*,10,Rainbow,1}
	ModifyGraph log(left)=1,log(vol)=1
	ModifyGraph mirror(left)=2,mirror(bottom)=2
	ModifyGraph font="Arial"
	ModifyGraph fSize=14
	ModifyGraph lblPos(left)=63,lblPos(vol)=61
	ModifyGraph lblLatPos(vol)=-10
	ModifyGraph freePos(vol)={0,bottom}
	ModifyGraph axisEnab(left)={0,0.48}
	ModifyGraph axisEnab(vol)={0.52,1}
	ModifyGraph dateInfo(bottom)={0,0,0}
	Label left "Diameter (nm)"
	Label bottom "Date and Time"
	Label vol "Volume dist."
	ColorScale/N=text0/A=RC/X=1.40/Y=-13.94/E image=smps_Nmatrix, heightPct=50
	ColorScale/C/N=text0 width=5, font="Arial", fsize=6
	ColorScale/N=text0_1/A=RC/X=3.37/Y=27.29/E image=smps_Vmatrix, heightPct=50
	ColorScale/C/N=text0_1 width=5, font="Arial", fsize=6
EndMacro

//************************************************************************************
Function PlotAverage_SMPS() 	:GraphMarquee
	string sdf = GetDataFolder(1)
	GetMarquee bottom
	setdatafolder StringByKey("YWAVEDF",imageinfo(S_marqueeWin,"smps_Nmatrix",0))
	
	variable i,j,k
	variable np,count
	wave date_time = $("smps_date_time")
	wave dia = $("smps_diameter")
	wave vmatrix = $("SMPS_Vmatrix")
	wave nmatrix = $("SMPS_Nmatrix")
	newdatafolder/o :avg_dist
//LIMIT on number of plots that can be made
//	variable limit_no = 500
//	for (i=0;i<limit_no;i+=1)
	variable avg_number
	i = 0
	do
		if (!waveexists($(":avg_dist:smps_avg_vol_dist"+num2str(i))))
			duplicate/o dia,$(":avg_dist:smps_avg_vol_dist"+num2str(i))
			duplicate/o dia,$(":avg_dist:smps_avg_num_dist"+num2str(i))
			wave avg_v = $(":avg_dist:smps_avg_vol_dist"+num2str(i))
			wave avg_n = $(":avg_dist:smps_avg_num_dist"+num2str(i))
			avg_number = i
			break
		endif
		i += 1
	while (1)
//	endfor

	avg_v = 0
	avg_n = 0
	np = numpnts(date_time)
	count = 0
	for (j=0;j<np;j+=1)
		if (date_time[j] > V_left && date_time[j] < V_right)
			if (numtype(vmatrix[j][0])==0)
				avg_v += vmatrix[j][p]
				avg_n += nmatrix[j][p]
				count+=1
			endif
		endif
	endfor
	avg_v /= count
	avg_n /= count
	dowindow/k Averaged_SMPS_plot
	Averaged_SMPS_plot_fxn(avg_n,avg_v,dia)
	
	//Determine the Scan #'s averaged
	variable a_pt = binarysearch(date_time,V_left)
	if (a_pt == -1) // Out of range, but placed before x_wave
		a_pt = 0
	elseif (a_pt == -2) //Out of range, but placed after x_wave
		a_pt = numpnts(date_time)
	endif
	if (date_time[a_pt] < V_left)
		a_pt +=1
	endif
	variable b_pt = binarysearch(date_time,V_right)
	if (b_pt == -1) // Out of range, but placed before x_wave
		b_pt = 0
	elseif (b_pt == -2) //Out of range, but placed after x_wave
		b_pt = numpnts(date_time)
	endif
	wave sample_no = $("smps_sample_no")
	variable start_scan = sample_no[a_pt]
	variable stop_scan = sample_no[b_pt]

	//Calculate and display statistics for selected averages
	wave dist_stats = calc_lognormal_statistics2(avg_n,dia)
	wave/z int_N = $(":avg_dist:int_N")
	wave/z mode_N = $(":avg_dist:mode_N")
	wave/z geom_N = $(":avg_dist:geom_N")
	wave/z sigma_N = $(":avg_dist:sigma_N")
	if (!waveexists(int_N))
		make/o/n=(avg_number+1) $(":avg_dist:int_N"), $(":avg_dist:mode_N"), $(":avg_dist:geom_N"), $(":avg_dist:sigma_N")
		wave int_N = $(":avg_dist:int_N")
		wave mode_N = $(":avg_dist:mode_N")
		wave geom_N = $(":avg_dist:geom_N")
		wave sigma_N = $(":avg_dist:sigma_N")	
	else
		redimension/n=(avg_number+1) int_N, mode_N, geom_N, sigma_N
	endif 
	int_N[avg_number] = dist_stats[0]
	mode_N[avg_number] = dist_stats[1]
	geom_N[avg_number] = dist_stats[2]
	sigma_N[avg_number] = dist_stats[3]

	wave dist_stats = calc_lognormal_statistics2(avg_v,dia)
	wave/z int_V = $(":avg_dist:int_V")
	wave/z mode_V = $(":avg_dist:mode_V")
	wave/z geom_V = $(":avg_dist:geom_V")
	wave/z sigma_V = $(":avg_dist:sigma_V")
	if (!waveexists(int_V))
		make/o/n=(avg_number+1) $(":avg_dist:int_V"), $(":avg_dist:mode_V"), $(":avg_dist:geom_V"), $(":avg_dist:sigma_V")
		wave int_V = $(":avg_dist:int_V")
		wave mode_V = $(":avg_dist:mode_V")
		wave geom_V = $(":avg_dist:geom_V")
		wave sigma_V = $(":avg_dist:sigma_V")	
	else
		redimension/n=(avg_number+1) int_V, mode_V, geom_V, sigma_V
	endif 
	int_V[avg_number] = dist_stats[0]
	mode_V[avg_number] = dist_stats[1]
	geom_V[avg_number] = dist_stats[2]
	sigma_V[avg_number] = dist_stats[3]

	string Stats_txt = "Average Dist:  Scan # " + num2str(start_scan) + " - " + num2str(stop_scan)+"\r"
	Stats_txt += "N_Int = " + num2str(int_N[avg_number]) + "\r"
	Stats_txt += "N_Mode = " + num2str(mode_N[avg_number]) + "\r"
	Stats_txt += "N_GeoMean = " + num2str(geom_N[avg_number]) + "\r"
	Stats_txt += "N_GeoStdev = " + num2str(sigma_N[avg_number])+"\r\r"
	Stats_txt += "V_Int = " + num2str(int_V[avg_number]) + "\r"
	Stats_txt += "V_Mode = " + num2str(mode_V[avg_number]) + "\r"
	Stats_txt += "V_GeoMean = " + num2str(geom_V[avg_number]) + "\r"
	Stats_txt += "V_GeoStdev = " + num2str(sigma_V[avg_number])
	TextBox/C/N=Nmode/A=MC Stats_txt

	killwaves dist_stats
	setdatafolder $sdf
End

//************************************************************************************
Function Statistics4SelectedMode_SMPS() 	:GraphMarquee
	string sdf = GetDataFolder(1)
	GetMarquee bottom
	//Marquee2Cursor() fxn----------------start
	string traces_on_graph = TraceNameList(S_marqueeWin, ";", 1+4)
	variable num_traces_on_graph = itemsinlist(traces_on_graph)
	variable index_of_trace
	SVAR /z Marquee2cursor_ChoosenTrace = root:Marquee2cursor_ChoosenTrace
	if (SVAR_Exists(Marquee2cursor_ChoosenTrace) && FindListItem(Marquee2cursor_ChoosenTrace,traces_on_graph)>=0 )
		traces_on_graph = RemoveFromList(Marquee2cursor_ChoosenTrace,traces_on_graph)
		traces_on_graph = Marquee2cursor_ChoosenTrace + ";" + traces_on_graph
	endif
	Prompt index_of_trace, "Choose Trace", popup, traces_on_graph
	DoPrompt "Choose Trace from Graph", index_of_trace	
	String /G root:Marquee2cursor_ChoosenTrace = Stringfromlist((index_of_trace-1),traces_on_graph)
	SVAR Marquee2cursor_ChoosenTrace = root:Marquee2cursor_ChoosenTrace
	wave y_wave = TraceNameToWaveRef(S_marqueeWin,Marquee2cursor_ChoosenTrace)
	wave x_wave = XWaveRefFromTrace(S_marqueeWin,nameofwave(y_wave))
	variable cursor_a_pt = binarysearch(x_wave,V_left)
	if (cursor_a_pt == -1) // Out of range, but placed before x_wave
		cursor_a_pt = 0
	elseif (cursor_a_pt == -2) //Out of range, but placed after x_wave
		cursor_a_pt = numpnts(x_wave)
	endif
	if (x_wave[cursor_a_pt] < V_left)
		cursor_a_pt +=1
	endif
	variable cursor_b_pt = binarysearch(x_wave,V_right)
	if (cursor_b_pt == -1) // Out of range, but placed before x_wave
		cursor_b_pt = 0
	elseif (cursor_b_pt == -2) //Out of range, but placed after x_wave
		cursor_b_pt = numpnts(x_wave)
	endif
//	print "Top YWave Name =", nameofwave(top_y_wave)
//	print "Top XWave Name =", nameofwave(top_x_wave)
//	print "V_left = ", V_left
//	print "V_right =", V_right
//	print "Cursor A Point =", cursor_a_pt
//	print "Cursor B Point =", cursor_b_pt
	Cursor /A=1/S=0/W=$S_marqueeWin A, $nameofwave(y_wave), cursor_a_pt
	Cursor /A=1/S=0/W=$S_marqueeWin B, $nameofwave(y_wave), cursor_b_pt
	//Marquee2Cursor() fxn----------------end	

	duplicate/FREE/R=[cursor_a_pt,cursor_b_pt] y_wave, y_dist
	duplicate/FREE/R=[cursor_a_pt,cursor_b_pt] x_wave, x_dist
	//Calculate and display statistics for selected averages
	wave dist_stats = calc_lognormal_statistics2(y_dist,x_dist)
	string SubStats_txt = "Mode Statistics: " + num2str(x_wave[cursor_a_pt]) + " - " + num2str(x_wave[cursor_b_pt]) + " nm\r"
	SubStats_txt += "INT = " + num2str(dist_stats[0]) + "\r"
	SubStats_txt += "MODE = " + num2str(dist_stats[1]) + "\r"
	SubStats_txt += "GEOMEAN = " + num2str(dist_stats[2]) + "\r"
	SubStats_txt += "GEOSTDEV = " + num2str(dist_stats[3])
	TextBox/A=MC SubStats_txt
	killwaves dist_stats
	setdatafolder $sdf
End

//************************************************************************************
function /WAVE calc_lognormal_statistics(dist,xwave)	//Calculate Lognormal Statistics for any given Distribution and X wave
	wave dist	
	wave xwave	
	//verify that dist and xwave have same dimensions
	if (numpnts(dist) != numpnts(xwave))
		Print "Problem - distribution and xwave have different dimensions.  Aborting."
	endif
	variable int=0, mode=0, geom=0, sigma=0, dlogX=0
	make/o/n=(numpnts(xwave)+1)/FREE xwave_im
	xwave_im = xwave[p]-(xwave[p+1] - xwave[p])/2
	xwave_im[numpnts(xwave_im)-2] = xwave[numpnts(xwave)-1] - (xwave[numpnts(xwave)-1] - xwave[numpnts(xwave)-2])/2
	xwave_im[numpnts(xwave_im)-1] = xwave[numpnts(xwave)-1] + (xwave[numpnts(xwave)-1] - xwave[numpnts(xwave)-2])/2
	wavestats/q dist
	mode = xwave[V_maxloc]
	variable i,j
	for (i=0;i<numpnts(dist);i+=1)
		dlogX = log(xwave_im[i+1]) - log(xwave_im[i])
		int += dist[i]*dlogX
		geom += dist[i]*dlogX*log(xwave[i])
	endfor
	geom = 10^(geom/int)
	for (j=0;j<numpnts(dist);j+=1)
		sigma += dist[j]*dlogX*(log(xwave[j])-log(geom))^2
	endfor	
	sigma = 10^((sigma/(int-1))^(1/2))
//	print "int = ", int
//	print "mode = ", mode
//	print "geom = ", geom
//	print "sigma = ", sigma
	make/o/n=4/FREE dist_stats
	dist_stats[0] = int		//integrated value of distribution
	dist_stats[1] = mode	//mode of distribution
	dist_stats[2] = geom	//geometric mean diameter (= Count Median Diameter for lognormal distributions)
	dist_stats[3] = sigma	//geometric standard deviation
	return dist_stats
end

//************************************************************************************
function /WAVE calc_lognormal_statistics2(dist,xwave)	//Calculate Lognormal Statistics for any given Distribution and X wave
	wave dist	//dX/dlogD, where X = N, S, V, etc.
	wave xwave
	//verify that dist and xwave have same dimensions
	if (numpnts(dist) != numpnts(xwave))
		Print "Problem - distribution and xwave have different dimensions.  Aborting."
	endif
	variable int=NaN,mode=NaN, geom=NaN, geos=NaN
	make/o/n=(numpnts(xwave)+1)/FREE xwave_im
	xwave_im = xwave[p]-(xwave[p+1] - xwave[p])/2
	xwave_im[numpnts(xwave_im)-2] = xwave[numpnts(xwave)-1] - (xwave[numpnts(xwave)-1] - xwave[numpnts(xwave)-2])/2
	xwave_im[numpnts(xwave_im)-1] = xwave[numpnts(xwave)-1] + (xwave[numpnts(xwave)-1] - xwave[numpnts(xwave)-2])/2
	make/d/n=(numpnts(xwave))/FREE lnDp,dlnDp, dist_x_dlnDp, dist_x_dlnDp_x_lnDp, dist_x_dlnDp_x_lnDp_m_lnGeom
	
	wavestats/q/M=1 dist
	if (V_npnts>0)
		mode = xwave[V_maxloc]
		lnDp = ln(xwave[p])
		dlnDp = ln(xwave_im[p+1]) - ln(xwave_im[p])
		dist_x_dlnDp = dist[p] * dlnDp[p] / 2.30258509299405		//NOTE:  dist uses LOG not LNb
		wavestats/q/M=1 dist_x_dlnDp
		int = V_sum
		dist_x_dlnDp_x_lnDp = dist_x_dlnDp[p] * lnDp[p]
		wavestats/q/M=1 dist_x_dlnDp_x_lnDp
		geom = exp(V_sum/int)			
		dist_x_dlnDp_x_lnDp_m_lnGeom = dist_x_dlnDp[p]*(lnDp[p]-ln(geom))^2
		wavestats/q/M=1 dist_x_dlnDp_x_lnDp_m_lnGeom
		geos = exp(sqrt(V_sum/int))	//Function for 1/N formulation
//		geos = exp(sqrt(V_sum/(int-1)))	//Function for 1/(N-1) formulation
	endif
//	print "int = ", int
//	print "mode = ", mode
//	print "geom = ", geom
//	print "sigma = ", sigma
	make/o/n=4 dist_stats
	dist_stats[0] = int		//integrated value of distribution
	dist_stats[1] = mode	//mode of distribution
	dist_stats[2] = geom	//geometric mean diameter (= Count Median Diameter for lognormal distributions)
	dist_stats[3] = geos	//geometric standard deviation
	return dist_stats
end

//Calculate Statistical information on PTOF image (dM/dlnDp vs Dp vs time)
//NOTE, originally wrote this in log()/10^() space, but IGOR had a mathematical problem somewhere, works in ln()/exp() space...
function calc_PTOFimage_stats_ln(PTOFimage, diameter_im) 

	wave PTOFimage	//dimensions = (np,nq)
//	wave date_time	//dimensions = np+1
	wave diameter_im		//dimensions = nq+1
	variable i,j
	
	variable np=dimsize(PTOFimage,0)
	variable nq=dimsize(PTOFimage,1)
	
	//Create diameter, lnDp, and dlnDp waves
	make/d/n=(nq)/FREE Dp, lnDp, dlnDp
	Dp = diameter_im[p]+(diameter_im[p+1]-diameter_im[p])/2
	Dp = exp(ln(diameter_im[p])+(ln(diameter_im[p+1])-ln(diameter_im[p]))/2)
	lnDp = ln(Dp[p])
	dlnDp=ln(diameter_im[p+1])-ln(diameter_im[p])
			
	//Make Calculation Waves = INT, MODE, GEOM, SIGMA
	make/d/n=(np)/FREE int=0, mode=NaN, geom=0, geos=0
	make/d/n=(nq)/FREE dist, dist_x_dlnDp, dist_x_dlnDp_x_lnDp, dist_x_dlnDp_x_lnDp_m_lnGeom

//	//GNU GSL Scientific Library Function - Function for 1/(N-1) formulation
//	make/d/n=(nq)/FREE dist_x_dlnDp_2
//	variable v0,v1
	
//Integrate IMAGE and calculate mode and geometric mean diameters
	for (i=0;i<np;i+=1)
		dist = PTOFimage[i][p]
		wavestats/q/M=1 dist
		if (V_maxloc>0)
			mode[i] = Dp[V_maxloc]
		else
			mode[i] = NaN
		endif
		if (V_npnts==0)
			int[i] = NaN
			geom[i] = NaN
			geos[i] = NaN
		else
			//Two calculation methods
	
			//FastOp and Wavestats
			dist_x_dlnDp = dist[p] * dlnDp[p]
			wavestats/q/M=1 dist_x_dlnDp
			int[i] = V_sum
			dist_x_dlnDp_x_lnDp = dist_x_dlnDp[p] * lnDp[p]
			wavestats/q/M=1 dist_x_dlnDp_x_lnDp
			geom[i] = exp(V_sum/int[i])			
			dist_x_dlnDp_x_lnDp_m_lnGeom = dist_x_dlnDp[p]*(lnDp[p]-ln(geom[i]))^2
			wavestats/q/M=1 dist_x_dlnDp_x_lnDp_m_lnGeom
			geos[i] = exp(sqrt(V_sum/int[i]))	//Function for 1/N formulation
			
//			//GNU GSL Scientific Library Function - Function for 1/(N-1) formulation
//			v0 = V_sum
//			dist_x_dlnDp_2 = (dist[p] * dlnDp[p])^2
//			wavestats/q/M=1 dist_x_dlnDp_2
//			v1 = V_sum
//			geos[i] = exp(sqrt((int[i]/(int[i]^2-v1))*v0))
	
//			//Standard using for-loops
//			for (j=0;j<nq;j+=1)
//				if ((numtype(dist[j])==0) && (numtype(dlnDp[j])==0))
//					int[i] += dist[j]*dlnDp[j]
//					if ((numtype(Dp[j])==0) && (Dp[j]>0))
//						geom[i] += dist[j]*dlnDp[j]*ln(Dp[j])
//					endif
//				endif
//			endfor
//			geom[i] = exp(geom[i]/int[i])
//			if (numtype(geom[i])==0)
//				for (j=0;j<nq;j+=1)
//					if ((numtype(dist[j])==0) && (numtype(dlnDp[j])==0))
//						if ((numtype(Dp[j])==0) && (Dp[j]>0))
//							geos[i] += dist[j]*dlnDp[j]*(ln(Dp[j])-ln(geom[i]))^2
//						endif
//					endif
//				endfor
//				geos[i] = exp(sqrt(geos[i]/int[i]))
//			endif
			
		endif	
	endfor
	
	//Set some pratical bounds
	variable max_diameter = 5000 	//nm
	variable min_diameter = 5 //nm
	mode = (mode[p] > max_diameter) || (mode[p] < min_diameter) ? (NaN):mode[p]
	geom = (geom[p] > max_diameter) || (geom[p] < min_diameter) ? (NaN):geom[p]
	variable max_sigma = 8
	geos = geos[p] > max_sigma ? (NaN):geos[p]
	
	//Make Output Wave
	make/o/n=(np,4)/d PTOFimage_stats
	SetDimLabel 1,0,Integration,PTOFimage_stats
	SetDimLabel 1,1,ModeDiameter,PTOFimage_stats
	SetDimLabel 1,2,GeometricMeanDiameter,PTOFimage_stats
	SetDimLabel 1,3,GeometricStandardDeviation,PTOFimage_stats
	PTOFimage_stats[][0] = int[p]
	PTOFimage_stats[][1] = mode[p]
	PTOFimage_stats[][2] = geom[p]
	PTOFimage_stats[][3] = geos[p]

	//Display output in Table
	edit /k=1 PTOFimage_stats
	modifytable horizontalIndex=2
	
	//Display output in Graph
	
end

//************************************************************************************
function Averaged_SMPS_plot_fxn(wave0,wave1,wave2)
	wave wave0,wave1,wave2
	Display /W=(340,200,720,430) wave0 vs wave2
	appendtograph/r wave1 vs wave2
	ModifyGraph lSize=2
	ModifyGraph grid(bottom)=1
	ModifyGraph log(bottom)=1
	ModifyGraph font="Arial"
	ModifyGraph fSize=14
	ModifyGraph standoff=0
	Label left "dN/dlogDp"
	Label right "dV/dlogDp"
	Label bottom "Dmobility (nm)"
	ModifyGraph axRGB(left)=(65280,0,0),axRGB(right)=(0,0,52224);DelayUpdate
	ModifyGraph tlblRGB(left)=(65280,0,0),tlblRGB(right)=(0,0,52224);DelayUpdate
	ModifyGraph alblRGB(left)=(65280,0,0),alblRGB(right)=(0,0,52224)
	execute("ModifyGraph rgb("+nameofwave(wave1)+")=(0,0,52224)")
end

//************************************************************************************
function fit_and_int_smps_1modeLognormal(low_size, high_size)
	variable low_size, high_size
	if (!waveexists(smps_date_time))
		print "==================================="
		print "No SMPS data found in current data folder.  Aborting."
		print "==================================="
		setdatafolder root:
		abort
	endif		
	wave dt = $("smps_date_time")
	wave dia = $("smps_diameter")
	wave geomean = $("smps_geomean_dia")
	wave vmatrix = $("smps_Vmatrix")
	wave int_v = $("smps_integrated_V")
	
	duplicate/o int_v, $("smps_int_V_FIT")
	wave int_v_fit = $("smps_int_V_FIT")
	Make/O/T/N=4 w_constraints
	wave/t w_constraints = w_constraints
	w_constraints[0] = {"K2 > 100","K2 < 500","K3 > 0.3","K3 < 1.5"}
	make/o/n=(dimsize(vmatrix,0),4) $("smps_V_fit_coefs")
	wave fit_coefs = $("smps_V_fit_coefs")
	variable i,j
	variable scale = 25

	duplicate/o dia, $("smps_log_diameter")
	wave log_dia = smps_log_diameter
	log_dia = log(dia)
	
	make/o/n = (dimsize(vmatrix,0)) fit_int_V
	wave fit_int_V
	make/o/n=1000 fit_dia, fit_log_dia
	wave fit_dia, fit_log_dia
	fit_dia = low_size+(high_size-low_size)/1000*p
	fit_log_dia = log(fit_dia)
	
	for (i=0;i<dimsize(vmatrix,0);i+=1)
		if (numtype(vmatrix[i][0])==0)
			make/o/n=(dimsize(vmatrix,1)) source
			wave source
			source = vmatrix[i][p]
			wavestats/q source
			scale = V_maxloc
			make/o/n=4 w_coefs = {0,scale,geomean[i],1}
			wave w_coefs
			CurveFit/Q/H="1000"/N LogNormal kwCWave=w_coefs,  source /X=dia /D ///C=w_constraints ///D=newfit 
			wave/z fit_source
			int_v_fit[i] = areaxy(log_dia, fit_source, -inf, inf)

			make/o/n=1000 fit_source2
			fit_source2 = w_coefs[0]+w_coefs[1]*exp(-(ln(fit_dia[p]/w_coefs[2])/w_coefs[3])^2)
			fit_int_V[i] = areaxy(fit_log_dia, fit_source2, -inf,inf)

			fit_coefs[i][] = w_coefs[q]
		endif
	endfor
	
	Integration_results_plot()
	CurveFit /Q line  smps_int_V_FIT /X=smps_integrated_V /D
	variable slope1 = K1 
	CurveFit /Q line  fit_int_V /X=smps_integrated_V /D
	variable slope2 = K1 
	Legend/N=text0/J/A=MC/X=-15.70/Y=42.21 "\\s(smps_integrated_V) smps_integrated_V"
	AppendText/N=text0 "\\s(smps_int_V_FIT) smps_int_V_FIT (slope = "+num2str(slope1)+")"
	AppendText/N=text0 "\\s(fit_int_V) fit_int_V (slope = "+num2str(slope2)+")"
	ShowInfo
	
	killwaves w_coef, w_constraints, w_paramconfidenceinterval, w_sigma
	killwaves source, fit_source, fit_source2
end

//************************************************************************************
// Fitting function if you want to fit 2 modes to a size distribution
Function LogNormal_2modes(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a0+a1*exp(-(ln(x/a2)/a3)^2)+b1*exp(-(ln(x/b2)/b3)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = a2
	//CurveFitDialog/ w[3] = a3
	//CurveFitDialog/ w[4] = b1
	//CurveFitDialog/ w[5] = b2
	//CurveFitDialog/ w[6] = b3

	return w[0]+w[1]*exp(-(ln(x/w[2])/w[3])^2)+w[4]*exp(-(ln(x/w[5])/w[6])^2)
End

//************************************************************************************
// Fitting function if you want to fit 2 modes to a size distribution
Function LogNormal_b(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (a1/(sqrt(2*pi)*log(a3)))*exp(-0.5*(ln(x/a2)/log(a3))^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = a2


	return (w[0]/(sqrt(2*pi)*log(w[2])))*exp(-0.5*(log(x/w[1])/log(w[2]))^2)
	//w[0]+w[1]*exp(-(ln(x/w[2])/w[3])^2)+w[4]*exp(-(ln(x/w[5])/w[6])^2)
End

//************************************************************************************
// Fitting function if you want to fit 2 modes to a size distribution
Function LogNormal_2modes_b(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a0+(a1/(sqrt(2*pi)*log(a3)))*exp(-0.5*(ln(x/a2)/log(a3))^2)+(b1/(sqrt(2*pi)*log(b3)))*exp(-0.5*(ln(x/b2)/log(b3))^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = a2
	//CurveFitDialog/ w[3] = a3
	//CurveFitDialog/ w[4] = b1
	//CurveFitDialog/ w[5] = b2
	//CurveFitDialog/ w[6] = b3

	return w[0]+(w[1]/(sqrt(2*pi)*log(w[3])))*exp(-0.5*(log(x/w[2])/log(w[3]))^2)+(w[4]/(sqrt(2*pi)*log(w[6])))*exp(-0.5*(log(x/w[5])/log(w[6]))^2)
	//w[0]+w[1]*exp(-(ln(x/w[2])/w[3])^2)+w[4]*exp(-(ln(x/w[5])/w[6])^2)
End

//************************************************************************************
// Fitting function if you want to fit 2 modes to a size distribution
Function LogNormal_3modes_b(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a0+(a1/(sqrt(2*pi)*log(a3)))*exp(-0.5*(ln(x/a2)/log(a3))^2)+(b1/(sqrt(2*pi)*log(b3)))*exp(-0.5*(ln(x/b2)/log(b3))^2)+(c1/(sqrt(2*pi)*log(c3)))*exp(-0.5*(ln(x/c2)/log(c3))^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = a2
	//CurveFitDialog/ w[3] = a3
	//CurveFitDialog/ w[4] = b1
	//CurveFitDialog/ w[5] = b2
	//CurveFitDialog/ w[6] = b3

	return w[0]+(w[1]/(sqrt(2*pi)*log(w[3])))*exp(-0.5*(log(x/w[2])/log(w[3]))^2)+(w[4]/(sqrt(2*pi)*log(w[6])))*exp(-0.5*(log(x/w[5])/log(w[6]))^2)+(w[7]/(sqrt(2*pi)*log(w[9])))*exp(-0.5*(log(x/w[8])/log(w[9]))^2)
	//w[0]+w[1]*exp(-(ln(x/w[2])/w[3])^2)+w[4]*exp(-(ln(x/w[5])/w[6])^2)
End

//************************************************************************************
// Fit and integrate a time series of size distributions
function fit_and_int_smps_2modeLognormal(low_size, high_size)
	variable low_size, high_size
	if (!waveexists(smps_date_time))
		print "==================================="
		print "No SMPS data found in current data folder.  Aborting."
		print "==================================="
		setdatafolder root:
		abort
	endif		
	wave dt = $("smps_date_time")
	wave dia = $("smps_diameter")
	wave geomean = $("smps_geomean_dia")
	wave vmatrix = $("smps_Vmatrix")
	wave int_v = $("smps_integrated_V")
	
	duplicate/o int_v, $("smps_int_V_FIT")
	wave int_v_fit = $("smps_int_V_FIT")
	Make/O/T/N=7 w_constraints
	wave/t w_constraints = w_constraints
	w_constraints[0] = {"K2 > 70","K2 < 300","K3 > 0.3","K3 < 1.5","K5 > 200","K5 < 800","K6 > 0.3","K6 < 1.5"}	//"K0 >(-0.1)","K0 < 10",
	make/o/n=(dimsize(vmatrix,0),7) $("smps_V_fit_coefs")
	wave fit_coefs = $("smps_V_fit_coefs")
	variable i,j
	variable scale

	duplicate/o dia, $("smps_log_diameter")
	wave log_dia = smps_log_diameter
	log_dia = log(dia)
	
	make/o/n = (dimsize(vmatrix,0)) fit_int_V
	wave fit_int_V
	make/o/n=1000 fit_dia, fit_log_dia
	wave fit_dia, fit_log_dia
	fit_dia = low_size+(high_size-low_size)/1000*p
	fit_log_dia = log(fit_dia)
	
	for (i=0;i<dimsize(vmatrix,0);i+=1)
		if (numtype(vmatrix[i][0])==0)
			make/o/n=(dimsize(vmatrix,1)) source
			wave source
			source = vmatrix[i][p]
			wavestats/q source
			scale = V_maxloc
			make/o/n=7 w_coefs = {0,scale/2,100,0.6,scale,400,0.6}
			wave w_coefs
			FuncFit/Q/H="1000000"/N LogNormal_2modes kwCWave=w_coefs,  source /X=dia /D /C=w_constraints ///D=newfit 
			wave/z fit_source
			int_v_fit[i] = areaxy(log_dia, fit_source)

			make/o/n=(numpnts(fit_dia)) fit_source2
			fit_source2 = w_coefs[0]+w_coefs[1]*exp(-(ln(fit_dia[p]/w_coefs[2])/w_coefs[3])^2)+w_coefs[4]*exp(-(ln(fit_dia[p]/w_coefs[5])/w_coefs[6])^2)
			fit_int_V[i] = areaxy(fit_log_dia, fit_source2, low_size, high_size)

			fit_coefs[i][] = w_coefs[q]
		endif
	endfor
	
	Integration_results_plot()
	CurveFit /Q line  smps_int_V_FIT /X=smps_integrated_V /D
	variable slope1 = K1 
	CurveFit /Q line  fit_int_V /X=smps_integrated_V /D
	variable slope2 = K1 
	Legend/N=text0/J/A=MC/X=-15.70/Y=42.21 "\\s(smps_integrated_V) smps_integrated_V"
	AppendText/N=text0 "\\s(smps_int_V_FIT) smps_int_V_FIT (slope = "+num2str(slope1)+")"
	AppendText/N=text0 "\\s(fit_int_V) fit_int_V (slope = "+num2str(slope2)+")"
	ShowInfo
	
	killwaves w_coef, w_constraints, w_paramconfidenceinterval, w_sigma
	killwaves source, fit_source, fit_source2
end

//************************************************************************************
Function Integration_results_plot() //: Graph
	Display /W=(75.75,330.5,905.25,704) smps_int_V_FIT,fit_int_V vs smps_integrated_V as "integration_results_plot"
	AppendToGraph/L=int/B=int_dt smps_integrated_V,smps_int_V_FIT, fit_int_V vs smps_date_time 
	ModifyGraph mode(smps_integrated_V)=4,mode(smps_int_V_FIT#1)=4,mode(fit_int_V#1)=4,mode(smps_int_V_FIT)=3
	ModifyGraph mode(fit_int_V)=3
	ModifyGraph marker(smps_integrated_V)=19,marker(smps_int_V_FIT)=19,marker(fit_int_V)=16
	ModifyGraph marker(smps_int_V_FIT#1)=19,marker(fit_int_V#1)=16
	ModifyGraph rgb(smps_int_V_FIT)=(0,0,0),rgb(fit_int_V)=(0,0,52224),rgb(smps_int_V_FIT#1)=(0,0,0)
	ModifyGraph rgb(fit_int_V#1)=(0,0,52224)
	ModifyGraph font="Arial"
	ModifyGraph fSize=14
	ModifyGraph standoff=0
	ModifyGraph lblPos(left)=63,lblPos(bottom)=54,lblPos(int)=60,lblPos(int_dt)=52
	ModifyGraph lblLatPos(left)=-9,lblLatPos(bottom)=-2,lblLatPos(int)=-12,lblLatPos(int_dt)=3
	ModifyGraph freePos(int)={0,int_dt}
	ModifyGraph freePos(int_dt)={0,int}
	ModifyGraph axisEnab(bottom)={0,0.45}
	ModifyGraph axisEnab(int_dt)={0.55,1}
	ModifyGraph dateInfo(bottom)={0,0,0}
	Label left "Integrated FIT Volume (um3/cm3) "
	Label bottom "Integrated SMPS Volume (um3/cm3) "
	Label int_dt "date and time"
	Label int "Integrated Volume (um3/cm3) "
//	Legend/N=text0/J/A=MC/X=-15.70/Y=42.21 "\\s(smps_integrated_V) smps_integrated_V\r\\s(smps_int_V_FIT) smps_int_V_FIT\r\\s(fit_int_V) fit_int_V"
End

function plot_dists1()
	if (strlen(CsrInfo(A)) <= 0)
		print "==================================="
		print "Cursor not found.  Aborting."
		print "==================================="
		abort
	endif		
	variable row = pcsr(A)
	print "row =", row
//	setdatafolder root:smps:
	wave vmatrix = smps_Vmatrix
	wave dia = smps_diameter
	wave fit_dia = fit_dia
	make/o/n= (dimsize(vmatrix,1)) dist,dist_fit
	make/o/n= (numpnts(fit_dia)) dist_fit2
	dist = vmatrix[row][p]
	wave coefs = smps_V_fit_coefs
	dist_fit = coefs[row][0]+coefs[row][1]*exp(-(ln(dia[p]/coefs[row][2])/coefs[row][3])^2)
	dist_fit2 = coefs[row][0]+coefs[row][1]*exp(-(ln(fit_dia[p]/coefs[row][2])/coefs[row][3])^2)
	
	dowindow/k Dists_plot
	display dist,dist_fit vs dia as "Dists_plot"
	dowindow/c Dists_plot
	appendtograph dist_fit2 vs fit_dia
	ModifyGraph rgb(dist_fit)=(0,0,0),rgb(dist_fit2)=(0,0,52224)
	ModifyGraph log(bottom)=1	
end

function plot_dists2()
	if (strlen(CsrInfo(A)) <= 0)
		print "==================================="
		print "Cursor not found.  Aborting."
		print "==================================="
		abort
	endif		
	variable row = pcsr(A)
	print "row =", row
//	setdatafolder root:smps:
	wave vmatrix = smps_Vmatrix
	wave dia = smps_diameter
	wave fit_dia = fit_dia
	make/o/n= (dimsize(vmatrix,1)) dist,dist_fit
	make/o/n= (numpnts(fit_dia)) dist_fit2
	dist = vmatrix[row][p]
	wave coefs = smps_V_fit_coefs
	dist_fit = coefs[row][0]+coefs[row][1]*exp(-(ln(dia[p]/coefs[row][2])/coefs[row][3])^2)+coefs[row][4]*exp(-(ln(dia[p]/coefs[row][5])/coefs[row][6])^2)
	dist_fit2 = coefs[row][0]+coefs[row][1]*exp(-(ln(fit_dia[p]/coefs[row][2])/coefs[row][3])^2)+coefs[row][4]*exp(-(ln(fit_dia[p]/coefs[row][5])/coefs[row][6])^2)

	dowindow/k Dists_plot
	display dist,dist_fit vs dia as "Dists_plot"
	dowindow/c Dists_plot
	appendtograph dist_fit2 vs fit_dia
	ModifyGraph rgb(dist_fit)=(0,0,0),rgb(dist_fit2)=(0,0,52224)
	ModifyGraph log(bottom)=1	
end

//************************************************************************************
function kill_list_of_waves_smps(wave_list)
	string wave_list
	variable i,j,k,ni
	ni = itemsinlist(wave_list)
	for (i=0;i<ni;i+=1)
		wave target = $(GetWavesDataFolder($stringfromlist(i,wave_list),2))
		killwaves target
	endfor
end

//************************************************************************************
function median_smps(waven)
	wave waven	
	variable result
	duplicate waven, tempMedianWave
	Sort tempMedianWave, tempMedianWave
	wavestats/q tempMedianWave
	redimension/n=(V_npnts) tempMedianWave
	SetScale/P x 0,1,tempMedianWave
	result = tempMedianWave((numpnts(tempMedianWave)-1)/2)
	KillWaves tempMedianWave
	return result
End

//************************************************************************************
function SplitIntoScans()
//&	setdatafolder root:smps
	string cdf = getdatafolder(1)	//Get Current Data Folder!

//TSI SMPS accounts for CPC flow dilution by comparing the input flowrate with the nominal instrument flow rate
	
//Concatenate all files together
	//find and concatenate all waves
	variable i,j
	variable ndf = countobjects(":",4)	//directories with files
	string datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("avg_dist") datafolder name from list
	if (stringmatch(datafolder_list,"*avg_dist*"))
		datafolder_list = RemoveFromList("avg_dist",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	
		
	//Lists of waves to be concatenated
	string cat_wave_list = "date_time_list;sample_no_list;time_diff_list;"
	string cat_wavename_list = "smps_date_time;smps_sample_no;smps_sample_time;"
	string/g date_time_list = ""
	string/g sample_no_list = ""
	string/g time_diff_list = ""
	//Lists of matrices to be concatenated
	string matrix_list = ""	
	
	for (i=0;i<ndf;i+=1)
		datafolder_str=stringfromlist(i,datafolder_list,",")	//note the string is separated by commas...	//GetIndexedObjName(":",4,j)
//		datafolder_str = GetIndexedObjName(":",4,i)

	//Check to see if "x" has been added to datafolder due to start of file name being a number...
		if (stringmatch(datafolder_str,"x*")==1)
			file_name_str = datafolder_str[1,strlen(datafolder_str)]
		else
			file_name_str = datafolder_str
		endif
//		if (stringmatch(datafolder_str,"avg_dist"))
//			ndf -= 1
//			continue
//		endif
		
	//Assume that all data uses same diameter wave, if NOT, then complain and quit...
//&		wave/z smps_diameter = $("root:smps:smps_diameter")
		wave/z smps_diameter = $(cdf+"smps_diameter")
		if (!waveexists(smps_diameter))		//(i==0)
//&			duplicate/o $("root:smps:"+datafolder_str+":dia_"+file_name_str), $("root:smps:smps_diameter")
//&			duplicate/o $("root:smps:"+datafolder_str+":dia_"+file_name_str+"_im"), $("root:smps:smps_diameter_im")
			duplicate/o $(cdf+datafolder_str+":dia_"+file_name_str), $("smps_diameter")
			duplicate/o $(cdf+datafolder_str+":dia_"+file_name_str+"_im"), $("smps_diameter_im")
		else
			//check if the diameter waves have changed between the data!!!!
//&			wave smps_diameter = $("root:smps:smps_diameter")
			wave smps_diameter = $(cdf+"smps_diameter")
//&			wave source_dia = $("root:smps:"+datafolder_str+":dia_"+file_name_str)
			wave source_dia = $(cdf+datafolder_str+":dia_"+file_name_str)
			if (!equalWaves(source_dia, smps_diameter, 1))
				print "Diameter waves are NOT the same!  Must process separately.  Quiting."
				abort
			endif
		endif

	//Generate lists of waves to be concatenated
//&		date_time_list += "root:smps:"+datafolder_str+":dt_"+file_name_str+";"
//&		sample_no_list += "root:smps:"+datafolder_str+":sn_"+file_name_str+";"
//&		time_diff_list += "root:smps:"+datafolder_str+":tdiff_"+file_name_str+";"
		date_time_list += cdf+datafolder_str+":dt_"+file_name_str+";"
		sample_no_list += cdf+datafolder_str+":sn_"+file_name_str+";"
		time_diff_list += cdf+datafolder_str+":tdiff_"+file_name_str+";"

	//Generate lists of matrices to be concatenated
//&		matrix_list += "root:smps:"+datafolder_str+":sdm_"+file_name_str+";"		
		matrix_list += cdf+datafolder_str+":sdm_"+file_name_str+";"		
		
	// Create folder to hold individual scans and split into scans
		string new_folder_str
		string scan_str
		string matrix_str
		matrix_str = cdf+datafolder_str+":sdm_"+file_name_str
		wave/z matrix_wv = $matrix_str
		if(!datafolderexists(cdf+datafolder_str+":scans"))
			new_folder_str = cdf+datafolder_str+":scans"
			newdatafolder $new_folder_str
		endif		
		variable nscans = dimsize($(matrix_str),0)
		for(j=0;j<=nscans;j+=1)
			scan_str = cdf+datafolder_str+":scans:Scan"+num2istr(j+1)
			make/o/d/n=(dimsize($matrix_str,1)) $Scan_str
			wave/z scan_wv = $scan_str
			scan_wv = matrix_wv[j][p]
		endfor
		duplicate/o $(cdf+datafolder_str+":dia_"+file_name_str), $(cdf+datafolder_str+":scans:smps_diameter")
	endfor
End

//************************************************************************************
Function Fit_Bimodal_byScan(file_name_str)
	string file_name_str

	variable i, j
	string cdf = getdatafolder(1)
	wave/z dia = $("smps_diameter")
	string scan_str
	
	string matrix_str = "::sdm_"+file_name_str
	variable nscans = dimsize($(matrix_str),0)
	variable dp_guess
	variable Np_guess
	variable width_guess = 0.25
	
	make/o/d/n=(numpnts(dia)) source
	wave source
	make/o/d/n=(nscans,7) fit_coefs = nan
	make/o/d/n=(nscans,numpnts(dia)) fit_curves = nan
	Make/O/T/N=7 w_constraints
	wave/t w_constraints = w_constraints
	w_constraints[0] = {"K2 > 70","K2 < 300","K3 > 0.3","K3 < 1.5","K5 > 200","K5 < 800","K6 > 0.3","K6 < 1.5"}	//"K0 >(-0.1)","K0 < 10",
	make/o/n=7 w_coefs

		for (i=0;i<nscans;i+=1)
		
			scan_str = "Scan"+num2istr(i+1)
			wave/z scan_wv = $scan_str
			source = scan_wv
			source = numtype(source)==2 ? 0 : source
			wavestats/q source
			dp_guess = dia[V_maxloc]
			Np_guess = source[V_maxloc]
			w_coefs = {0,Np_guess/10,dp_guess*0.8,width_guess,Np_guess,dp_guess,width_guess}
			wave w_coefs
			w_constraints[0] = {"K1 > 10","K1 < " + num2str(Np_guess*1.1),"K2 > 5","K2 < "+num2str(dp_guess*1),"K3 > 0.1","K3 < 0.5","K4 > 10","K4 < "+num2str(Np_guess*1.1),"K5 > 5","K5 < "+num2str(dp_guess*1.3),"K6 > 0.1","K6 < 0.5"}
			FuncFit/Q/H="1000000"/N LogNormal_2modes kwCWave=w_coefs,  source /X=dia /D /C=w_constraints ///D=newfit 
			wave fit_source
			
			if(width_guess==w_coefs[6])
				w_coefs = nan
				fit_source = nan
			endif
			fit_coefs[i][] = w_coefs[q]
			fit_curves[i][] = w_coefs[0]+w_coefs[1]*exp(-(ln(dia[q]/w_coefs[2])/w_coefs[3])^2)+w_coefs[4]*exp(-(ln(dia[q]/w_coefs[5])/w_coefs[6])^2)
		endfor
end

//************************************************************************************
// Function to merge dXdlogDp from multiple SMPS files into one long time series
// Currently not general b/c must set folders to be used (must update "fldrs_str" wave)
// Typically run this after loading the AIM data
// Merged t-series stored in "merged" folder
Function MergeSMPS_tseries(datatype)
	string datatype // e.g. dNdlogDp, dN, dNdDp, etc.

	string df_SMPS = "root:SMPS:"
	setdatafolder $df_SMPS
	string df_SMPS_Merged = "root:SMPS:Merged"
	if(datafolderexists(df_SMPS_Merged)==0)
		newdatafolder $df_SMPS_Merged
	endif
	
	variable ndf = countobjects(":",4)	//directories with files
	string datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("Merged") datafolder name from list
	if (stringmatch(datafolder_list,"*Merged*"))
		datafolder_list = RemoveFromList("Merged",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	

	setdatafolder df_SMPS_merged

	make/o/n=(ndf) wv_npnts = nan
	
	variable V_npntstot 
	variable V_dia
	
	string cdf_str
	string wvstr_dia1 = "dia_"
	string wvstr_time = "dt_"
	string current_str
	
	variable i, j
		
	make/o/d/n=(V_dia) dia_SMPS = nan
	wave dia_SMPS
	
	make/o/d/n=(V_dia+1) dia_SMPS_im = nan
	wave dia_SMPS_im
	
	string dt_list = "", dt_im_list="", sdm_list=""
	
	for(i=0;i<ndf;i+=1)
		cdf_str = df_SMPS + stringfromlist(i,datafolder_list,",")//[i]
		setdatafolder $cdf_str
		
		current_str = wvstr_dia1 + stringfromlist(i,datafolder_list,",")
		wave current = $current_str
		if(i==0) 
			duplicate/o current, $(df_SMPS+"merged:dia_SMPS")//dia_SMPS = current
		endif
		
		current_str = "dia_" + stringfromlist(i,datafolder_list,",") + "_im"
		wave current = $current_str
		if(i==0)
			duplicate/o current, $(df_SMPS+"merged:dia_SMPS_im")//dia_SMPS_im = current
		endif
		
		dt_list += df_SMPS+stringfromlist(i,datafolder_list,",")+":"
		dt_list += wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		
		if(i<ndf-1)
			dt_im_list += df_SMPS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		else
			dt_im_list += df_SMPS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += wvstr_time+stringfromlist(i,datafolder_list,",")+"_im;"
		endif
		
		sdm_list += df_SMPS+stringfromlist(i,datafolder_list,",")+":"
		sdm_list += "sdm_"+stringfromlist(i,datafolder_list,",")+";"

	endfor
	
	setdatafolder df_SMPS_Merged	
	concatenate/o dt_list, dt_SMPS
	setscale/p d 0,1, "dat", dt_SMPS
	
	concatenate/o dt_im_list, dt_SMPS_im
	setscale/p d 0,1, "dat", dt_SMPS_im
	
	concatenate/o/np=0 sdm_list, sdm_SMPS
	duplicate/o sdm_SMPS, $(datatype+"_SMPS")
	wave sdm_merged = $(datatype+"_SMPS")
	killwaves sdm_SMPS
	
	sdm_merged = numtype(sdm_merged)==2 ? 0 : sdm_merged
	
	killwaves/z wv_npnts
	
	setdatafolder df_SMPS_Merged	
End

//************************************************************************************
// Function to merge multiple SMPS files into one long time series
// Currently not general b/c must set folders to be used (must update "fldrs_str" wave)
Function MergeAPS_tseries(datatype)
	string datatype // e.g. dNdlogDp, dN, dNdDp, etc.

	string df_APS = "root:APS:"
	setdatafolder $df_APS
	string df_APS_Merged = "root:APS:Merged"
	if(datafolderexists(df_APS_Merged)==0)
		newdatafolder $df_APS_Merged
	endif
	
	variable ndf = countobjects(":",4)	//directories with files
	string datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("Merged") datafolder name from list
	if (stringmatch(datafolder_list,"*Merged*"))
		datafolder_list = RemoveFromList("Merged",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	

	setdatafolder df_APS_merged

	make/o/n=(ndf) wv_npnts = nan
	
	variable V_npntstot 
	variable V_dia
	
	string cdf_str
	string wvstr_dia1 = "dia_"
	string wvstr_time = "dt_"
	string current_str
	
	variable i, j
	
	make/o/d/n=(V_dia) dia_APS = nan
	wave dia_APS
	
	make/o/d/n=(V_dia+1) dia_APS_im = nan
	wave dia_APS_im
	
	string dt_list = "", dt_im_list="", sdm_list=""
	
	for(i=0;i<ndf;i+=1)
		cdf_str = df_APS + stringfromlist(i,datafolder_list,",")//[i]
		setdatafolder $cdf_str
		
		current_str = wvstr_dia1 + stringfromlist(i,datafolder_list,",")
		wave current = $current_str
		if(i==0) 
			duplicate/o current, $(df_APS+"merged:dia_APS")//dia_APS = current
		endif
		
		current_str = "dia_" + stringfromlist(i,datafolder_list,",") + "_im"
		wave current = $current_str
		if(i==0)
			duplicate/o current, $(df_APS+"merged:dia_APS_im")//dia_APS_im = current
		endif
		
		dt_list += df_APS+stringfromlist(i,datafolder_list,",")+":"
		dt_list += wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		
		if(i<ndf-1)
			dt_im_list += df_APS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		else
			dt_im_list += df_APS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += wvstr_time+stringfromlist(i,datafolder_list,",")+"_im;"
		endif
		
		sdm_list += df_APS+stringfromlist(i,datafolder_list,",")+":"
		sdm_list += "sdm_"+stringfromlist(i,datafolder_list,",")+";"

	endfor
	
	setdatafolder df_APS_Merged	
	concatenate/o dt_list, dt_APS
	setscale/p d 0,1, "dat", dt_APS
	
	concatenate/o dt_im_list, dt_APS_im
	setscale/p d 0,1, "dat", dt_APS_im
	
	concatenate/o/np=0 sdm_list, sdm_APS
	duplicate/o sdm_APS, $(datatype+"_APS")
	wave sdm_merged = $(datatype+"_APS")
	killwaves sdm_APS
	
	sdm_merged = numtype(sdm_merged)==2 ? 0 : sdm_merged
	
	killwaves/z wv_npnts
	setdatafolder df_APS_Merged	
End

//**************************************************************************************************
// assumes particular timing
// Ambient = starts on zeros
// TD = starts on 5's
Function SplitSMPS_byTDstate()
	
//	setdatafolder root:smps:Merged
	string dt_smps_str="dt_smps"
	string dndlogdp_str="dNdlogDp_smps"
	variable tref_ambient=0
	variable tref_denuded=5
	variable tref_treatspecial1=25
	variable tref_treatspecial2=55
	prompt dt_smps_str, "Select time wave", popup WaveList("*",";","")
	prompt dndlogdp_str, "Select number distribution wave",popup WaveList("*",";","")
	prompt tref_ambient, "Starting minute of ambient data"
	prompt tref_denuded, "Starting minute of TD data"
	prompt tref_treatspecial1, "Minute where only 1/2 a scan was completed"
	prompt tref_treatspecial2, "Minute where only 1/2 a scan was completed"
	DoPrompt "Select waves", dt_smps_str, dndlogdp_str,tref_ambient,tref_denuded, tref_treatspecial1, tref_treatspecial2
	
	wave dt_SMPS = $dt_smps_str
	wave dndlogdp_SMPS = $dndlogdp_str
	
	if(numpnts(dt_smps)!=dimsize(dndlogdp_smps,0))
		Abort "Waves must have same length. Aborting."
	endif
	
	variable V_time = numpnts(dt_SMPS)
	make/o/d/n=(V_time) maskTD = nan, maskSpecial=1
	variable i
	string dt_string, dt, dtt
	
	for(i=0;i<V_time;i+=1)
		dt_string = secs2time(dt_smps[i],2)
		dt = dt_string[4]//stringfromlist(1,dt_string,":")
		dtt = dt_string[3,4]
		if(str2num(dt)==tref_ambient)
			maskTD[i] = 0
		elseif(str2num(dt)==tref_denuded)
			maskTD[i] = 1
		endif
		if(str2num(dtt)==tref_treatspecial1 || str2num(dtt)==tref_treatspecial2)
			maskSpecial[i] = 2
		else
			maskSpecial[i] = 1
		endif
	endfor
	
	string TD = dndlogdp_str + "_TD"
	string RT = dndlogdp_str + "_RT"
	duplicate/o dndlogdp_SMPS $TD, $RT
	wave dndlogdp_TD = $TD
	wave dndlogdp_RT = $RT
	dNdlogDp_TD[][] = maskTD[p]==1 ? dndlogdp_SMPS[p][q] : nan
	dNdlogDp_RT[][] = maskTD[p]==0 ? dndlogdp_SMPS[p][q] : nan
	dNdlogDp_TD[][] *= maskSpecial[p]
	
	killwaves/z maskTD, maskSpecial
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Macro to allow for averaging of matrices (i.e. size distributions)
// adapted from AverageOnlyXSecs from 2010GeneralMacros.ipf
// cdcappa@ucdavis.edu, 27 January 2013
Macro AveragSizeDistOnlyXSecs(timeWaveStr, wave2avgStr, waveAvgStr, XSec)
	String timeWaveStr, wave2avgStr, waveAvgStr= "my_avg"
	Variable XSec=60
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup WaveList("*",";","")
	prompt wave2avgStr, "Matrix Data to be averaged?",popup WaveList("*",";","")
	prompt waveAvgStr,  "Name of the averaged wave?"
	prompt XSec, "Seconds to average over?"

	Silent(1)
	
	if (numpnts($timeWaveStr)!=dimsize($wave2avgStr,0))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragOnlyXSecs"
	endif	
		
	if (strlen(waveAvgStr)==0)
		Abort "You need to insert values into the name prompts - Aborting from AveragOnlyXSecs"
	endif
		
	variable startTime =$timeWaveStr[0] - mod($timeWaveStr[0], XSec)
	variable stopTime = $timeWaveStr[numpnts($timeWaveStr)-1] -  mod($timeWaveStr[numpnts($timeWaveStr)-1], XSec ) + XSec
	
	Make/o/d/n=((stopTime-startTime)/XSec) $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	SetScale d 0,0,"dat", $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec")
	$("start_"+num2str(XSec)+"sec") = startTime+p*XSec
	$("stop_"+num2str(XSec)+"sec") = startTime+(p+1)*XSec - 1
	
	// for image plotting
	Make/o/d/n=((stopTime-startTime)/XSec+1) $("start_"+num2str(XSec)+"sec_im"), $("stop_"+num2str(XSec)+"sec_im"); DelayUpdate
	SetScale d 0,0,"dat", $("start_"+num2str(XSec)+"sec_im"), $("stop_"+num2str(XSec)+"sec_im")
	$("start_"+num2str(XSec)+"sec_im") = startTime+p*XSec; DelayUpdate
	$("stop_"+num2str(XSec)+"sec_im") = startTime+(p+1)*XSec - 1; DelayUpdate
	
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=((stopTime-startTime)/XSec,dimsize($wave2avgStr,1)) $waveAvgStr; DelayUpdate
	
	Make/o/d/n=(dimsize($wave2avgStr,0)) wave2Avg; DelayUpdate
	Make/o/d/n=((stopTime-startTime)/XSec) waveAvgCol; DelayUpdate
	variable i=0
	variable nd = dimsize($wave2avgStr,1)
	
	do
		wave2Avg = $wave2AvgStr[p][i]
		DoAveragOnlyUsingStartStop($timeWaveStr, wave2Avg, waveAvgCol, $("start_"+num2str(XSec)+"sec"), $("stop_"+num2str(XSec)+"sec"))
		$waveAvgStr[][i] = waveAvgCol[p];DelayUpdate
		i+=1
	while(i<nd)
	
	killwaves/z wave2avg, waveAvgCol
	DoUpdate
	
End	

// **************************************************************************************************************************************
// Macro to allow for averaging of matrices (i.e. size distributions) to selected start/stop periods
// adapted from AverageOnlyXSecs from 2010GeneralMacros.ipf
// cdcappa@ucdavis.edu, 27 January 2013
Macro AveragSizeDistUsingStartStop(timeWaveStr, wave2avgStr, waveAvgStr,startWaveStr, stopWaveStr)
	String timeWaveStr, wave2avgStr, waveAvgStr="my_avg",startWaveStr, stopWaveStr
	prompt timeWaveStr, "Timewave of the data to be averaged?",popup WaveList("*",";","")
	prompt wave2avgStr, "Data to be averaged?",popup WaveList("*",";","")
	prompt waveAvgStr, "Name to give the averaged wave?"
	prompt startWaveStr, "Which wave has the start times?",popup WaveList("*",";","")
	prompt stopWaveStr, "Which wave has the stop times?",popup WaveList("*",";","")

	if (numpnts($timeWaveStr)!=dimsize($wave2avgStr,0))
		Abort " The first two waves must have the same number of points."+timeWaveStr + " "+wave2avgStr+" - Aborting from AveragOnlyXSecs"
	endif	
	
	if (numpnts($startWaveStr)!=numpnts($stopWaveStr))
		Abort " The start and stop waves must have the same number of points."+ startWaveStr+" "+stopWaveStr + " - Aborting from AveragOnlyUsingStartStop"
	endif	
	
	if (strlen(waveAvgStr)==0) 
		Abort "You need to insert values into the name prompts - Aborting from AveragOnlyUsingStartStop"
	endif
	
	silent(1)
	
	// for image plotting
	Make/o/d/n=(numpnts($startWaveStr)+1) $(startWaveStr+"_im")
	SetScale d 0,0,"dat", $(startWaveStr+"_im")
	$(startWaveStr+"_im") = MakeTraceForImage(startWaveStr) 
		
	waveAvgStr = CleanUpName(waveAvgStr, 0)
	Make/o/d/n=(Numpnts($startWaveStr), dimsize($wave2avgStr,1)) $waveAvgStr
	
	Make/o/d/n=(dimsize($wave2avgStr,0)) wave2Avg
	Make/o/d/n=(numpnts($startWaveStr)) waveAvgCol
	variable i=0
	variable nd = dimsize($wave2avgStr,1)

	do
		wave2Avg = $wave2AvgStr[p][i]
		DoAveragOnlyUsingStartStop($timeWaveStr, wave2Avg, waveAvgCol, $startWaveStr, $stopWaveStr)
		$waveAvgStr[][i] = waveAvgCol[p]
		i+=1
	while(i<nd)
	
	killwaves/z wave2avg, waveAvgCol
 
 End
	

//*************************************************************************************
// Interpolate a given size distribution time-series onto a different time-base
// Can be used to "fill in" missing periods
// Requires having "2010GeneralMacros.ipf" loaded
Function Interp_SizeDist_tseries(st_dNdlogDp,st_time,st_dp,st_time_interp)
	string st_dNdlogDp // matrix with size distribution (time is columns, rows are size)
	string st_time // time wave of original distribution
	string st_dp // diameter wave
	string st_time_interp // start time for averaging
	
	string cdf = getdatafolder(2) // current data folder
	string st_dNdlogDp_orig, st_time_orig, st_dp_orig
	string st_dNdlogDp_interp
	
	st_dNdlogDp_orig = cdf+st_dNdlogDp
	st_time_orig = cdf+st_time
	st_dp_orig = cdf+st_dp
	st_dNdlogDp_interp = cdf + st_dNdlogDp + "_interp" 
	
	wave time_interp = $st_time_interp
	wave dpwave = $st_dp
	wave dndlogdp = $st_dNdlogDp_orig
	wave timew = $st_time_orig
	
	// check to make sure you are in the folder you think and that you have the waves you need
	if(waveexists(dndlogdp)==0 || waveexists(dpwave)==0 || waveexists(timew)==0 || waveexists(time_interp)==0)// || waveexists(stop_time)==0)
		print "Missing some important waves...are you sure you're in the right folder?"
		abort
	endif
	
	variable npnts_new = numpnts(time_interp) // number of time steps in new size distribution matrix
	variable npnts_old = numpnts(timew) // number of time steps in original size distribution matrix
	variable nsize = numpnts(dpwave) // number of diameters
	
	make/o/d/n=(npnts_old) dN_single // a wave to hold a time series of dNdlogDp for each size (i.e. an array)
	make/o/d/n=(npnts_new,nsize) $st_dNdlogDp_interp // a new matrix to hold the final averaged dNdlogDp
	wave dndlogdp_interp = $st_dNdlogDp_interp
	duplicate/o/d dndlogdp dndlogdp_temp // duplicate original wave, so as to get rid of nans
	dndlogdp_temp = numtype(dndlogdp)==2 ? 0 : dndlogdp
	
	variable i
	
	for(i=0;i<nsize;i+=1)
		dn_single = dndlogdp_temp[p][i]
		// linearly interpolate waves onto desired time base
		Interpolate2/T=1/N=200/I=3/Y=dN_single_L/X=time_interp timew, dN_single
		wave dN_single_L
		dNdlogDp_interp[][i] = dN_single_L[p]

	endfor
	
	duplicate/o dndlogdp_interp dSdlogDp_interp, dVdlogDp_interp
	dsdlogdp_interp = pi*(dpwave[q]^2)*dndlogdp_interp
	dvdlogdp_interp = (pi/6)*(dpwave[q]^3)*dndlogdp_interp
	
	// make sure to have "extended" waves necessary for making image plots
	string st_time_interp_i = st_time_interp + "_i"
	make/o/d/n=(npnts_new+1) $st_time_interp_i
	wave time_interp_i = $st_time_interp_i
	time_interp_i = time_interp
	time_interp_i[npnts_new] = time_interp_i[npnts_new-1]+(time_interp_i[npnts_new-1]-time_interp_i[npnts_new-2])
	setscale d, 0, 1, "dat", time_interp_i
	
	string st_dp_i = st_dp + "_i"
//	if(waveexists($st_dp_i)==0)
		make/o/d/n=(nsize+1) $st_dp_i
		wave dp_i = $st_dp_i
		dp_i = dpwave
		dp_i[nsize] = dp_i[nsize-1] + (dp_i[nsize-1]-dp_i[nsize-2])
//	endif
	
	killwaves/z dn_single, dn_single_L my_avg, my_num, my_std, dNdlogDp_temp
	setdatafolder $cdf
	
End

//*************************************************************************************
// Interpolate a given size distribution time-series onto a different time-base
// Can be used to "fill in" missing periods
// Requires having "2010GeneralMacros.ipf" loaded
Function InterpSMPSforTD()//st_dNdlogDp,st_time,st_dp,st_time_interp)

	string st_dNdlogDp="dNdlogDp" // matrix with size distribution (time is columns, rows are size)
	
	prompt st_dndlogdp, "Select size distribution wave", popup WaveList("*",";","")
	DoPrompt "Select waves" st_dndlogdp
	
	string st_dNdlogDp_interp
	st_dNdlogDp_interp = st_dndlogdp + "_interp" 
	
	wave dndlogdp = $st_dndlogdp
	variable nt = dimsize(dndlogdp,0)
	variable nd = dimsize(dndlogdp,1)
	
	make/o/d/n=(nd) dn1,dn2,dn3//dN_single // a wave to hold a time series of dNdlogDp for each size (i.e. an array)
	make/o/d/n=(nt,nd) $st_dNdlogDp_interp // a new matrix to hold the final averaged dNdlogDp
	wave dndlogdp_interp = $st_dNdlogDp_interp
	dndlogdp_interp = dndlogdp
	
	variable i
	
	for(i=0;i<nt;i+=1)
		if(numtype(dndlogdp[i][0])==2)
			dn1 = dndlogdp[i-1][p]
			dn2 = dndlogdp[i+1][p]
			dn3 = (dn1+dn2)/2
			dndlogdp_interp[i][] = dn3[q]
		endif
	endfor
	
	killwaves/z dn1,dn2,dn3
	
End


//***************************************************************************************************************************************
// Graph size distribution time series
Function GraphSizeDist_tseries(dXdlogDp_st,tseries_st,dp_st,scalingfactor)
	string dXdlogDp_st
	string tseries_st
	string dp_st
	variable scalingfactor
	
	wave dndlogdp = $dXdlogDp_st
	wave tseries = $tseries_st
	wave dp = $dp_st
	
	// check to make sure you are in the folder you think and that you have the waves you need
	if(waveexists(dndlogdp)==0 || waveexists(dp)==0 || waveexists(tseries)==0)
		print "Missing some important waves...are you sure you're in the right folder?"
		abort
	endif

	Display;AppendImage $dXdlogDp_st vs {$tseries_st,$dp_st};DelayUpdate
	wavestats/q $dXdlogDp_st
	ModifyImage $dXdlogDp_st ctab= {0,scalingfactor*v_max,Rainbow,0}
	ModifyGraph log(left)=1;DelayUpdate
	Label left "Particle Diameter (nm)";DelayUpdate
	Label bottom "Date and Time"
	TextBox/C/N=text1 dXdlogDp_st

End

//***************************************************************************************************************************************
// Graph size distribution time series
Function GraphSizeDist_tseries_pnts(dXdlogDp_st,dp_st,start,stop,scalingfactor)
	string dXdlogDp_st
	string dp_st
	variable start
	variable stop
	variable scalingfactor
	
	wave dndlogdp = $dXdlogDp_st
	wave dp = $dp_st
	
	// check to make sure you are in the folder you think and that you have the waves you need
	if(waveexists(dndlogdp)==0 || waveexists(dp)==0)
		print "Missing some important waves...are you sure you're in the right folder?"
		abort
	endif

	Display
	variable i
	for(i=start;i<=stop;i+=1)
		AppendToGraph dNdlogDp[i][] vs $dp_st;DelayUpdate
		ModifyGraph log(bottom)=1;DelayUpdate
		Label bottom "Particle Diameter (nm)";DelayUpdate
		Label left "Concentration"
	endfor

End

////*********************************************************************************************************************
//// convert size distributions between different forms
//Macro Convert_dNdlogDp(dNStr, DpStr, dSstr, dVStr)
//	String dNStr, DpStr, dSstr="dSdlogDp", dVstr="dVdlogDp"
//	prompt dNStr, "Name of the number-weighted size distribution time series to be averaged?",popup WaveList("*",";","")
//	prompt DpStr, "Name of the Particle diameter wave?",popup WaveList("*",";","")
//	prompt dSStr, "Name to give the surface area-weighted wave?"
//	prompt dVStr, "Name to give the volume-weighted wave?"
//
//	if (dimsize($dNStr,1)!=numpnts($DpStr))
//		Abort " The diameter does not match with the time series."+dNStr + " "+DpStr+ " - Aborting from Number2SandV_tseries"
//	endif	
//	
//	if (strlen(dVStr)==0 || strlen(dSstr)==0) 
//		Abort "You need to insert values into the name prompts - Aborting from Number2SandV_tseries"
//	endif
//		
//	dSstr = CleanUpName(dSstr, 0)
//	Make/o/n=(dimsize($dNstr,0), dimsize($dNstr,1)) $dSstr
//	$dSstr[][] = (4*pi)*($dpstr[q]/2)^2*$dNstr[p][q]
//	
//	dVstr = CleanUpName(dVstr, 0)
//	Make/o/n=(dimsize($dNstr,0), dimsize($dNstr,1)) $dVstr
//	$dVstr[][] = (pi/6)*$dpstr[q]^3*$dNstr[p][q]
//	
// End
 
 //*********************************************************************************************************************
// convert size distributions between different forms
Function Convert_dNdlogDp(pr)//dNStr, DpStr, dSstr, dVStr)
	variable pr
	
	String dNStr="dNdlogDp", DpStr="Dp", dSstr="dSdlogDp", dVstr="dVdlogDp"
	prompt dNStr, "Name of the number-weighted size distribution time series to be converted?",popup WaveList("*",";","")
	prompt DpStr, "Name of the Particle diameter wave?",popup WaveList("*",";","")
	prompt dSStr, "Name to give the surface area-weighted wave?"
	prompt dVStr, "Name to give the volume-weighted wave?"
	if(pr==1)
		DoPrompt "Wave Selection", dNstr, DpStr, dSstr, dVstr
	endif

	wave dNdlogDp = $dNstr
	wave Dp = $DpStr
	
	if (dimsize(dNdlogDp,1)!=numpnts(Dp))
		Abort " The diameter does not match with the time series."+dNStr + " "+DpStr+ " - Aborting from Number2SandV_tseries"
	endif	
	
	if (strlen(dVStr)==0 || strlen(dSstr)==0) 
		Abort "You need to insert values into the name prompts - Aborting from Number2SandV_tseries"
	endif
		
	dSstr = CleanUpName(dSstr, 0)
	Make/o/d/n=(dimsize(dNdlogDp,0), dimsize(dNdlogDp,1)) $dSstr
	wave dSdlogDp = $dSstr
	dSdlogDp[][] = (4*pi)*(Dp[q]*1e-9/2)^2*dNdlogDp[p][q]*1e12
	
	dVstr = CleanUpName(dVstr, 0)
	Make/o/d/n=(dimsize(dNdlogDp,0), dimsize(dNdlogDp,1)) $dVstr
	wave dVdlogDp = $dVstr
	dVdlogDp[][] = (pi/6)*(Dp[q]*1e-9)^3*dNdlogDp[p][q]*1e18
	
 End

//******************************************************
 Function Number2Volume(dNdlogDp,Dp)
 	wave dndlogdp
 	wave dp
 	
 	make/o/n=(numpnts(dp)) dVdlogDp_temp
 	dVdlogdp_temp = (pi/6)*dp^3*dndlogdp
 	
// 	return dVdlogdp_temp
 	
 End
 	
//********************************************************************************************
Function ConvertDpaToDpm()
	variable density=2.0
	string APSstr, DpaStr
	variable type=0, DpaUnits=0
	string MobilityStr="dNdlogDpm", DpmStr="Dpm"	
	prompt density, "What density do you want to use? "
	prompt APSstr, "Select your dN or dNdlogDp time series", popup WaveList("*",";","")
	prompt type, "Is your APS wave dN or dNdlogDp?", popup, "dN;dNdlogDp"
	prompt DpaStr, "Select your aerodynamic diameter wave", popup WaveList("*",";","")
	prompt DpaUnits, "Is your diameter in microns or nm?",popup, "microns;nm"
	prompt MobilityStr, "What do you want to name your converted size distribution?"
	prompt dpmStr, "What do you want to name your mobility diameter wave?"
	DoPrompt "Set Density", density, APSstr, type, DpaStr, DpaUnits, MobilityStr, dpmStr
	if (V_Flag)
		return 0									// user canceled
	endif
	
	if(dimsize($APSstr,1)!=numpnts($DpaStr))
		Abort "Size distribution and particle diameter waves must have the same dimensions. Aborting ConvertDpaToDpm."
	endif
	
	wave APSwave = $APSstr
	wave dpa = $dpastr
	variable nd = numpnts(dpa)
	variable nt = dimsize(APSwave,0)
	make/o/n=(nt,nd) $mobilitystr
	wave MobilityWave = $mobilityStr
	make/o/n=(nd) $dpmstr
	wave Dpm = $dpmstr
	make/o/n=(nd+1) $(dpmstr+"_im")
	wave Dpm_im = $(dpmstr+"_im")
	
	make/o/n=(nd)/FREE APS_single
	make/o/n=(nd)/FREE logDpa=log(Dpa), logDpm
	make/o/n=(nd)/FREE slip_a, slip_m=1
	make/o/n=(nd)/FREE dpa2
	
	if(dpaunits==1)
		dpa2 = dpa*1000
	elseif(dpaunits==2)
		dpa2 = dpa
	endif
	
	// Slip Correction
	slip_a =1+2*(0.065/(dpa2/1000))*(1.257+0.4*EXP(-0.55*(dpa2/1000)/0.065))
	// loop through to find new dpm
	variable i=0
	do
		dpm = dpa2*sqrt((slip_a/slip_m)*(1/density))
		slip_m =1+2*(0.065/(dpm/1000))*(1.257+0.4*EXP(-0.55*(dpm/1000)/0.065))
	i+=1
	while(i<100)
	
	note/K dpm, "Dpa converted to Dpm assuming a density of " + num2str(density) + "; units = nm"
	note/K dpm_im, "Dpa converted to Dpm assuming a density of " + num2str(density) + "; units = nm"
	note/K mobilitywave, "Dpa converted to Dpm assuming a density of " + num2str(density)
	
	dpm_im = dpm
	dpm_im[numpnts(dpm_im)] = dpm[numpnts(dpm)] + (dpm[numpnts(dpm)-1]-dpm[numpnts(dpm)-2])
	
	// Assume that the user wants output as dNdlogDp, not dN
	logDpm = log(dpm)
	Differentiate logdpm/D=dlogdpm
	
	if(type==1) 										// dN = dN
		mobilitywave = APSwave / dlogDpm[q]
	elseif(type==2) 									// convert dNdlogDp to dN
		Differentiate logDpa/D=dlogDpa
		mobilitywave = APSwave*(dlogDpa[q]/dlogDpm[q])
	endif
	
	killwaves/z dlogDpa, dlogDpm, logDpa, logDpm, slip_a, slip_m, dpa2
	
end

//********************************************************************************************
Function ConvertDpaToDpmSingle(dpa,density)
	variable dpa // units = nm
	variable density

	variable slip_a, slip_m=1
	variable dpa2 = dpa
	variable dpm
		
	// Slip Correction
	slip_a =1+2*(0.065/(dpa2/1000))*(1.257+0.4*EXP(-0.55*(dpa2/1000)/0.065))
	// loop through to find new dpm
	variable i=0
	do
		dpm = dpa2*sqrt((slip_a/slip_m)*(1/density))
		slip_m =1+2*(0.065/(dpm/1000))*(1.257+0.4*EXP(-0.55*(dpm/1000)/0.065))
	i+=1
	while(i<100)
	
	return dpm
	
end

//********************************************************************************************
// First set data folder to the folder containing aerodynamic diameter waves
// This function will assume that your density wave has been created and is in this folder
// requires as input APS data as dN (not dNdlogDp)
Function ConvertDpaToDpm_DensityWave()
	string APSstr="dN_APS_avg", DpaStr="dia_APS"
	variable type=0, DpaUnits=0
	string MobilityStr="dNdlogDpm_APS", DpmStr="Dpm_APS"	
	string DensityStr="Density_fromAMS"
	prompt densityStr, "Select wave containing Density time series", popup WaveList("*",";","")
	prompt APSstr, "Select your dN or dNdlogDp time series", popup WaveList("*",";","")
	prompt type, "Is your APS wave dN or dNdlogDp?", popup, "dN;dNdlogDp"
	prompt DpaStr, "Select your aerodynamic diameter wave", popup WaveList("*",";","")
	prompt DpaUnits, "Is your diameter in microns or nm?",popup, "microns;nm"
	prompt MobilityStr, "What do you want to name your converted size distribution?"
	prompt dpmStr, "What do you want to name your mobility diameter wave?"
	DoPrompt "Set Density", densityStr, APSstr, type, DpaStr, DpaUnits, MobilityStr, dpmStr
	if (V_Flag)
		return 0									// user canceled
	endif
	
	if(dimsize($APSstr,1)!=numpnts($DpaStr))
		Abort "Size distribution and particle diameter waves must have the same dimensions. Aborting ConvertDpaToDpm."
	endif
	
	wave APSwave = $APSstr
	wave Density = $DensityStr
	wave dpa = $dpastr
	variable nd = numpnts(dpa)
	variable nt = dimsize(APSwave,0)
	make/o/d/n=(nt,nd) $mobilitystr
	wave MobilityWave = $mobilityStr
	make/o/d/n=(nd) $dpmstr
	wave Dpm = $dpmstr
	make/o/d/n=(nd+1) $(dpmstr+"_im")
	wave Dpm_im = $(dpmstr+"_im")
	
	make/o/n=(nd)/FREE APS_single
	make/o/d/n=(nd)/FREE Dpm_single, logDpm_single, dN_single
	make/o/n=(nd)/FREE logDpa=log(Dpa), logDpm
	make/o/n=(nd)/FREE slip_a, slip_m=1
	make/o/n=(nd)/FREE dpa2
	make/o/d/n=(nt,nd) Dpm_matrix=nan, dN_matrix=nan
	
	if(dpaunits==1)
		dpa2 = dpa*1000
	elseif(dpaunits==2)
		dpa2 = dpa
	endif
	
	variable i,j
	// Get matrix of Dpm values based on density wave and Dpa
	for(j=0;j<nt;j+=1)
		if(numtype(Density[i])!=2)
			// Slip Correction
			slip_a =1+2*(0.065/(dpa2/1000))*(1.257+0.4*EXP(-0.55*(dpa2/1000)/0.065))
			// loop through to find new dpm
			i=0
			do
				dpm_single = dpa2*sqrt((slip_a/slip_m)*(1/density[j]))
				slip_m =1+2*(0.065/(dpm_single/1000))*(1.257+0.4*EXP(-0.55*(dpm_single/1000)/0.065))
			i+=1
			while(i<20)
			Dpm_matrix[j][] = dpm_single[q]
		endif
	endfor

	// Dpm = average from individual Dpm's
	make/o/n=(nt) Dpm_time
	for(i=0;i<nd;i+=1)
		Dpm_time = Dpm_matrix[p][i]
		wavestats/q Dpm_time
		Dpm[i] = V_avg
	endfor
	
	// Interpolate all of the dN values with their individual Dpm's to the average Dpm
	for(i=0;i<nt;i+=1)
		dpm_single = Dpm_matrix[i][p]
		dN_single = APSwave[i][p]
		wavestats/q dN_single
		if(V_numnans==0)
			Interpolate2/T=1/I=3/Y=dN_L/X=dpm_single Dpm, dN_single
			dN_matrix[i][] = dN_L[q]
		endif
	endfor
	
	note/K dpm, "Dpa converted to Dpm using Density Wave; units = nm"
	note/K dpm_im, "Dpa converted to Dpm using Density Wave; units = nm"
	note/K mobilitywave, "Dpa converted to Dpm using Density Wave; units = nm"
	
	dpm_im = dpm
	dpm_im[numpnts(dpm_im)] = dpm[numpnts(dpm)] + (dpm[numpnts(dpm)-1]-dpm[numpnts(dpm)-2])
	
	// Assume that the user wants output as dNdlogDp, not dN
	logDpm = log(dpm)
	Differentiate logdpm/D=dlogdpm
	
	if(type==1) 										// dN = dN
		mobilitywave = dN_matrix / dlogDpm[q]
	elseif(type==2) 									// convert dNdlogDpa to dNdlogDpm
		Differentiate logDpa/D=dlogDpa
		mobilitywave = dN_matrix*(dlogDpa[q]/dlogDpm[q])
	endif
	
	killwaves/z dlogDpa, dlogDpm, logDpa, logDpm, slip_a, slip_m, dpa2, dpm_single, logdpm_single
	killwaves/z dN_single, dN_matrix
	
end

 //********************************************************************************************
 // Must set to appropriate folder first
Macro ConvertDpaToDpmMacro(density,APSstr,DpaStr,type,DpaUnits,MobilityStr,DpmStr)
//Function ConvertDpaToDpm()
	variable density=2.0
	string APSstr, DpaStr
	variable type=0, DpaUnits=0
	string MobilityStr="dNdlogDpm", DpmStr="Dpm"	
	prompt density, "What density do you want to use? "
	prompt APSstr, "Select your dN or dNdlogDp time series", popup WaveList("*",";","")
	prompt type, "Is your APS wave dN or dNdlogDp?", popup, "dN;dNdlogDp"
	prompt DpaStr, "Select your aerodynamic diameter wave", popup WaveList("*",";","")
	prompt DpaUnits, "Is your diameter in microns or nm?",popup, "microns;nm"
	prompt MobilityStr, "What do you want to name your converted size distribution?"
	prompt dpmStr, "What do you want to name your mobility diameter wave?"
	
	if(dimsize($APSstr,1)!=numpnts($DpaStr))
		print "Size distribution and particle diameter waves must have the same dimensions. Aborting ConvertDpaToDpm."
		abort
	endif
	
	Silent(1)
	
	variable nd = numpnts($dpastr)
	variable nt = dimsize($APSstr,0)
	make/o/n=(nt,nd) $mobilitystr
	make/o/n=(nd) $dpmstr
	make/o/n=(nd+1) $(dpmstr+"_im")
	
	make/o/n=(nd) logDpa=log($DpaStr), logDpm
	make/o/n=(nd) slip_a, slip_m=1
	make/o/n=(nd) dpa2
	
	if(dpaunits==1)
		dpa2 = $(dpastr)*1000
	else
		dpa2 = $(dpastr)
	endif
	
	// Convert Dpa --> Dpm
	slip_a =1+2*(0.065/(dpa2/1000))*(1.257+0.4*EXP(-0.55*(dpa2/1000)/0.065))					// Cunningham slip correction factor, aerodynamic diameter
	variable i=0
	do
		$(dpmstr) = dpa2*sqrt((slip_a/slip_m)*(1/density))										// mobility diameter from aerodynamic diameter
		slip_m =1+2*(0.065/($(dpmstr)/1000))*(1.257+0.4*EXP(-0.55*($(dpmstr)/1000)/0.065))		// Cunningham slip correction factor, mobility diameter
	i+=1
	while(i<100)
	
	note/K $dpmstr, "Dpa converted to Dpm assuming a density of " + num2str(density) + "; units = nm"
	note/K $(dpmstr+"_im"), "Dpa converted to Dpm assuming a density of " + num2str(density) + "; units = nm"
	note/K $mobilitystr, "Dpa converted to Dpm assuming a density of " + num2str(density)
	
	$(dpmstr+"_im") = $dpmstr
	$(dpmstr+"_im")[numpnts($(dpmstr+"_im"))] = $dpmstr[numpnts($dpmstr)] + ($dpmstr[numpnts($dpmstr)-1]-$dpmstr[numpnts($dpmstr)-2])
	
	// Assume that the user wants output as dNdlogDp, not dN
	logDpm = log($dpmstr)
	Differentiate logdpm/D=dlogdpm						// width of size bin, mobility diameter
	
	if(type==1) 										
		$mobilitystr = $APSstr / dlogDpm[q]				// convert dN to dNdlogDpm
	else									
		Differentiate logDpa/D=dlogDpa						// width of size bin, aerodynamic diameter
		$mobilitystr = $APSstr*(dlogDpa[q]/dlogDpm[q])		// convert dNdlogDpa to dNdlogDpm
	endif
	
	killwaves/z dlogDpa, dlogDpm, logDpa, logDpm, slip_a, slip_m, dpa2
	
end

 //********************************************************************************************
Function DpmShapeCorrection(Dpm,SCF)
	// convert mobility diameters to shape-corrected mass diameters
	variable Dpm // = mobility diameter, in nm
	variable SCF // = shape correction factor; AmmSulf = 1.02; NaCl = 1.08
	variable slip_mob
	variable slip_mass
	variable Dpmass = Dpm
	variable MFP = 68 // nm, mean free path
	
	// Convert Dpmobility --> Dpmass
	slip_mob =1+2*(MFP/Dpm)*(1.142+0.558*EXP(-0.55*Dpm/MFP))				// Cunningham slip correction factor, aerodynamic diameter
	slip_mass = slip_mob
	variable i=0
	do
		slip_mass =1+2*(MFP/dpmass)*(1.142+0.558*EXP(-0.55*dpmass/MFP))		// Cunningham slip correction factor, mobility diameter
		Dpmass = dpm*(slip_mass/slip_mob)*(1/SCF)										// mobility diameter from aerodynamic diameter
	i+=1
	while(i<100)
	
	return(DpMass)
	
end

//*******************************************************************************************************************
Function MakeTraceForImage(wvstr)
	string wvstr
	wave wv = $wvstr
	
	make/o/d/n=(numpnts(wv)+1) wv_im
	make/o/d/n=(numpnts(wv))/FREE dwv = wv
	differentiate wv /D=dwv
	wv_im = wv
	wv_im[numpnts(wv_im)] = wv_im[numpnts(wv_im)] + dwv[numpnts(wv)]
	
//	return wv_im
End

//*****************************************************************************************************************
// Merge SMPS and APS traces at some particular diameter
// assume diameter is in nm
// assumes data have been averaged to same time base with equally spaced times, but that the don't necessarily start/stop at the same time
Function MergeSMPSwithAPS(dpmerge,upperlimit)
	variable dpmerge				// diameter at which to merge
	variable upperlimit 			// upper limit APS diameter to include

 	variable GetNames = 1
 	
 	string df_merge = "root:SizeDistMerged"
 	if(datafolderexists(df_merge)==0)
 		newdatafolder/s $df_merge
 	endif
 		setdatafolder $df_merge
 		
//	variable density = 1.5
	string df_smps_str = "root:smps:merged"
	string df_aps_str = "root:aps:merged"

	prompt dpmerge, "Diameter to merge?"
	prompt upperlimit, "Upper limit APS diameter?"
//	prompt density, "Assumed Density?"
	prompt df_smps_str, "SMPS data folder"
	prompt df_aps_str, "APS data folder"
 	DoPrompt "Merge Diameter", dpmerge, upperlimit, df_smps_str, df_aps_str
 	Print "SMPS and APS data merged using a split at " + num2str(dpmerge) + " nm"
 	Print " and an upper limit cutoff of " + num2str(upperlimit) + " nm"
			
	DFREF df_smps = $df_smps_str
	DFREF df_aps = $df_aps_str
	
	string dpm_smps_str = "dia_smps"
	string dndlogdp_smps_str = "dndlogdp_smps_avg_RT_interp"
	string starttime_smps_str = "starttime"
	string dpm_aps_str = "dpm_aps"
	string dndlogdp_aps_str = "dndlogdpm_aps_avg_scale"
	string starttime_aps_str = "starttime"
	
	if(GetNames==1)
		prompt dpm_smps_str, "Dpm SMPS wave"
		prompt dndlogdp_smps_str, "dNdlogDp SMPS wave"
		prompt starttime_smps_str, "SMPS Time Wave"
		prompt dpm_aps_str, "Dpm APS wave"
		prompt dndlogdp_aps_str, "dNdlogDpm APS wave"
		prompt starttime_aps_str, "APS Time Wave"
		DoPrompt "Wave Names", dpm_smps_str, dpm_aps_str, dndlogdp_smps_str, dndlogdp_aps_str, starttime_smps_str, starttime_aps_str
	endif
			
	wave dpm_smps = df_smps:$dpm_smps_str
	wave dndlogdp_smps = df_smps:$dndlogdp_smps_str
	wave starttime_smps = df_smps:$starttime_smps_str
	wave dpm_aps = df_aps:$dpm_aps_str
	wave dndlogdp_aps = df_aps:$dndlogdp_aps_str
	wave starttime_aps = df_aps:$starttime_aps_str
	string aps_note = note(dndlogdp_aps)	
	
	variable nd_smps = numpnts(dpm_smps)
	variable nd_aps = numpnts(dpm_aps)
	variable nt_smps = dimsize(dndlogdp_smps,0)
	variable nt_aps = dimsize(dndlogdp_aps,0)
	
		variable idex_smps
	if(dpmerge>dpm_smps[nd_smps])
		idex_smps = nd_smps-1
	else
		findlevel/Q/P dpm_smps, dpmerge
		idex_smps = floor(V_LevelX)
	endif
	findlevel/Q/P dpm_aps, dpmerge
	variable idex_aps = floor(V_levelX)
	wavestats/q dpm_aps
	variable idex_aps_UL
	if(upperlimit > V_max) // dpm_aps)
		idex_aps_UL = numpnts(dpm_aps)
	else
		findlevel/Q/P dpm_aps, upperlimit
		idex_aps_UL = floor(V_levelX) + 1
	endif
	
	variable dpcut_smps = dpm_smps[idex_smps]
	variable dpcut_aps = dpm_aps[idex_aps]
	if(dpcut_smps > dpcut_aps)
		idex_aps += 1								// APS diameter must be greater than SMPS diameter
	endif
	
	variable ndcut_smps = idex_smps
	variable ndcut_aps = nd_aps - idex_aps + 1 - (nd_aps - idex_aps_UL)
	variable nd_merge = ndcut_smps + ndcut_aps
	
	variable tstart_smps = starttime_smps[0]
	variable tstop_smps = starttime_smps[numpnts(starttime_smps)]
	variable tstart_aps = starttime_aps[0]
	variable tstop_aps = starttime_aps[numpnts(starttime_aps)]
	wavestats/q starttime_smps
	variable deltatime = (V_max-V_min)/(V_npnts-1)
	variable tstart
	variable tstop
	variable tindex_smps
	variable tindex_aps
	
	if(tstart_smps<tstart_aps)
		tstart = starttime_smps[0]
		tindex_smps = 0
		findvalue/z /V=(tstart) starttime_aps
		tindex_aps = V_value
	else
		tstart = starttime_aps[0]
		tindex_aps = 0
		findvalue/z /V=(tstart) starttime_smps
		tindex_smps = V_value
	endif
	
	if(tstop_smps<tstop_aps)
		tstop = starttime_smps[numpnts(starttime_smps)]
	else
		tstop = starttime_aps[numpnts(starttime_aps)]
	endif
	
	variable nt_merge = (tstop-tstart)/deltatime + 1					// number of time steps in merged wave
	
	make/o/d/n=(nt_merge) tseries									// time series
	setscale d 0,0, "dat", tseries
	tseries = tstart + p*deltatime
	make/o/d/n=(nt_merge+1) tseries_im							// time series for image plotting
	setscale d 0,0, "dat", tseries_im
	tseries_im = tstart + p*deltatime
	
	make/o/d/n=(nt_merge,nd_merge) dNdlogDp = nan			// dNdlogDp for merged size distribution
	note/K dNdlogDp, "p/cm^3"
	make/o/d/n=(nt_merge,nd_merge) dSdlogDp = nan			// dSdlogDp for merged size distribution
	note/K dSdlogDp, "um^2 per m^3"
	make/o/d/n=(nt_merge,nd_merge) dVdlogDp = nan			// dVdlogDp for merged size distribution
	note/K dVdlogDp, "um^3 per m^3"
//	make/o/d/n=(nt_merge,nd_merge) dMdlogDp = nan			// dMdlogDp for merged size distribution
//	note/K dMdlogDp, "ug per m^3; Assumes density = " + num2str(density)
	make/o/d/n=(nd_merge) Dp = nan							// diameter for merged size distribution
	note/K Dp, aps_note
	make/o/d/n=(nd_merge+1) Dp_im = nan						// diameter for image plot of merged size distribution
	note/K Dp_im, aps_note
	Dp[0,ndcut_smps] = Dpm_smps[x]
	Dp[ndcut_smps+1,] = Dpm_aps[x-ndcut_smps+idex_aps-1]
	Dp_im = MakeTraceForImage("Dp")
	
	make/o/n=(nd_smps) logDp_smps = log(Dpm_smps)
	make/o/n=(nd_aps) logDp_aps = log(Dpm_aps)
	differentiate logDp_smps /D=dlogDp_smps
	differentiate logDp_aps /D=dlogDp_aps
	
	dNdlogDp[tindex_smps,][0,ndcut_smps] = dNdlogDp_smps[p-tindex_smps][q]
	dndlogDp[tindex_aps,][ndcut_smps+1,] = dNdlogDp_aps[p-tindex_aps][q-ndcut_smps+idex_aps-1]
	
	dSdlogDp = 4*pi*((Dp[q]*1e-9/2)^2)*dNdlogDp[p][q]*1e12				// dSdlogDp in microns squared per m^3
	dVdlogDp = (pi/6)*((Dp[q]*1e-9)^3)*dNdlogDp[p][q]*1e18				// dVdlogDp in um^3 per m^3											// g/cm^3
//	dMdlogDp = 1e6 * dVdlogDp * density								// dMdlogDp in micrograms per m^3
	
	killwaves/z dlogdp_smps, dlogdp_aps, ddp, logDp_smps, logDp_aps

End

//**********************************************************************************
function Diurnal_SizeDist(timewave,datawave,newwavestr,normalize,[maskwave])
	wave timewave, datawave
	string newwavestr
	variable normalize
	wave maskwave
	
	variable nt=(numpnts(timewave))
	variable nd = dimsize(datawave,1)
	
	if(paramisdefault(maskwave))
		make/o/d/n=(nt)/FREE masktemp = 1
		wave maskwave = masktemp
	endif
	
	if(waveexists(timewave)==0 || waveexists(datawave)==0)
		Abort "One or more waves do not exist"
	elseif(numpnts(timewave) != nt)
		Abort "Data wave is not as long as time wave"
	endif
	
	// must have time wave in Julian Day
	if(timewave[0]>400) 
		
		variable YearVar, YearSecs, SecsPerDay=60*60*24
		make/o/d/n=(nt)/FREE twave
		string datestr = secs2date(timewave[0],-2)
		sscanf datestr, "%i", yearvar
		YearSecs = date2secs(yearvar,1,1)
		twave = (TimeWave - YearSecs)/SecsPerDay
	else
		make/o/d/n=(nt)/FREE twave
		twave = timewave
	endif
	
	variable nhours = 24
	variable i, j
	variable time_i, time_f
	
	make/o/d/n=(nt) dN
	make/o/d/n=(nhours,nd) $newwavestr
	wave newwave = $newwavestr

	string str = newwavestr + "_sem"
	make/o/d/n=(nhours,nd) $str
	wave newwave_sem = $str
	
	make/o/d/n=(nhours) Time_diurnal
	make/o/d/n=(nhours+1) Time_diurnal_im
	make/o/d/n=(nd)/FREE dNdlogDp_Single
	
	// Determine averages by hour of day
	for(j=0;j<nd;j+=1)
		dN = datawave[p][j]
		dN = maskwave==0 ? nan : dN
	
		for(i=0;i<nhours;i+=1)
			time_i = i/24
			time_f = (i+1)/24
		extract dN, diurnalwave, (twave[x] - floor(twave[x])) >= time_i && (twave[x] - floor(twave[x])) < time_f
		newwave[i][j] = meanwithnan("diurnalwave",0)
		wavestats/q diurnalwave
		newwave[i][j] = V_avg
		newwave_sem[i][j] = V_sem
		Time_diurnal[i] = (Time_i+time_f)*24/2
		endfor
	endfor
	
	// Normalize size distributions to 1
	if(normalize==1)
		for(i=0;i<nhours;i+=1)
			dNdlogDp_single = newwave[i][p]
			wavestats/q dNdlogDp_single
			newwave[i][] = dNdlogDp_single[q]/V_max
		endfor
	endif
	
	Time_diurnal_im = Time_diurnal
	Time_Diurnal_im[numpnts(Time_diurnal_im)] = Time_diurnal[numpnts(Time_diurnal)] + 1
	
	killwaves/z diurnalwave, dN
end

//***********************************************************************************
Function DisplayDiurnalSizeDist(dXdlogDp, dpwave,[appendit])
	wave dXdlogDp, dpwave
	variable appendit // 0 = new graph
	
	if(ParamIsDefault(appendit))
		appendit = 0
	endif
	variable nt = dimsize(dXdlogDp,0)
	variable i
	
	if(appendit==0)
		display
	endif
	for(i=0;i<nt;i+=1)
		appendtograph dXdlogDp[i][] vs dpwave
	endfor
end

//***********************************************************************************
Function AppendDiurnalSizeDist(dXdlogDp, dpwave)
	wave dXdlogDp, dpwave
	
	variable nt = dimsize(dXdlogDp,0)
	variable i
	
	for(i=0;i<nt;i+=1)
		appendtograph dXdlogDp[i][] vs dpwave
	endfor
end

//**********************************************************************************
function BinByWave_SizeDist(binwavestr,datawave,newwavestr,normalize, start, step, nBins,[logscale,maskwave])
	string binwavestr // name of wave to which you want to bin
	wave datawave // name of the size distribution wave
	string newwavestr  // name of the new wave
	variable normalize // do you want to normalize
	variable start // start value; actually this value minus 0.5*step size
	variable step // step size
	variable nBins // number of bins
	variable logscale // use log scale, in which case step size is log-based
	wave maskwave // 0 = exclude, 1 = include
	
	if(paramisdefault(logscale))
		logscale=0
	endif
	
	wave binwave = $binwavestr
	
	variable nt=(numpnts(binwave))
	variable nd = dimsize(datawave,1)
	
	if(waveexists(binwave)==0 || waveexists(datawave)==0)
		Abort "One or more waves do not exist"
	elseif(dimsize(datawave,0) != nt)
		Abort "Data wave is not as long as time wave"
	endif
	
	variable i, j
	
	make/o/d/n=(nt) dN
	make/o/d/n=(nBins,nd) $newwavestr
	wave newwave = $newwavestr
	setscale/p x, start, step, newwave
	note/K newwave "Binned to " + binwavestr
	string str = newwavestr + "_sem"
	make/o/d/n=(nBins,nd) $str
	wave newwave_sem = $str
	
	make/o/d/n=(nBins+1) BinWaveX
	BinWaveX = start + step*x - step/2
	if(logscale==1)
		BinWaveX = 10^(log(start) + step*x - step/2)
	endif
	
	make/o/d/n=(nd)/FREE dNdlogDp_Single
	
	// Determine averages by binning
	for(j=0;j<nd;j+=1)
		dN = datawave[p][j]
		if(paramisdefault(maskwave)==0)
			dN = maskwave==1 ? dN : nan
		endif
	
		for(i=0;i<nBins;i+=1)
			extract dN, tempwave, binwave[x] >= (BinWaveX[i]) && binwave[x] < (BinWaveX[i]+step)
			if(numpnts(tempwave)>0)
				newwave[i][j] = meanwithnan("tempwave",0)
				wavestats/q tempwave
				newwave[i][j] = V_avg
				newwave_sem[i][j] = V_sem
			else
				newwave[i][j] = nan
				newwave_sem[i][j] = nan
			endif
		endfor
	endfor
	
	// Normalize size distributions to 1
	if(normalize==1)
		for(i=0;i<nBins;i+=1)
			dNdlogDp_single = newwave[i][p]
			wavestats/q dNdlogDp_single
			newwave[i][] = dNdlogDp_single[q]/V_max
		endfor
	endif
		
	killwaves/z tempwave, dN//, BinWaveX
end

//***********************************************************************************
Function IntegrateSizeDist([CutDiameter])
	variable CutDiameter // set this value if you want to integrate to some diameter < max

	if(paramisdefault(CutDiameter))
		CutDiameter = -1
	endif
	
	 variable density = 1.5
	 string df_merge = getdatafolder(2)//"root:SizeDistMerged"
	 string dndlogdp_str = "dNdlogDp"
	 string dsdlogdp_str = "dSdlogDp"
	 string dVdlogdp_str = "dVdlogDp"
	 string dMdlogDp_str = "dMdlogDp"
	 string dp_str = "Dp"
	 
	 prompt density, "Density g/cm^3"
	 prompt df_merge, "Data folder"
	 prompt dndlogdp_str, "dNdlogDp", popup WaveList("*",";","")
	 prompt dsdlogdp_str, "dSdlogDp", popup WaveList("*",";","")
	 prompt dVdlogdp_str, "dVdlogDp", popup WaveList("*",";","")
	 prompt dp_str, "Dp", popup WaveList("*",";","")
	 DoPrompt "Data folder and waves", density, df_merge, dndlogdp_str, dsdlogdp_str, dvdlogdp_str, dp_str
	 
	 setdatafolder $df_merge
	 
	 wave dndlogdp = $dndlogdp_str
	 wave dsdlogdp = $dsdlogdp_str
	 wave dvdlogdp = $dvdlogdp_str
	 wave dp = $dp_str
	 variable nd = dimsize(dndlogdp,1)
	 variable nt = dimsize(dndlogdp,0)
	 variable CutIndex
	 
	 if(paramisdefault(CutDiameter))
		CutDiameter = dp[nd-1]
		CutIndex = nd-1
	else
		FindLevel/Q/P dp, CutDiameter
		if(V_Flag==0)
			CutIndex = ceil(V_LevelX)-1
			CutDiameter = dp[CutIndex]
		else
			abort "You'd better double check your cut diameter, fool."
		endif
	endif
	 
	 if(waveexists(dsdlogdp)==0 || waveexists(dVdlogDp)==0)
	 	Print "Need to convert dNdlogDp to dSdlogDp and dVdlogDp"
	 	convert_dndlogdp(1)
	 endif
	 
	 make/o/d/n=(nt) Number_Conc, SurfaceArea, Volume_Conc, Mass_conc
	 duplicate/o dp logDp
	 logDp = log(dp)
	 
	 integrate/DIM=1/METH=1 dNdlogDp /X=logDp /D=dXdlogDp_int
	 Number_conc = dxdlogdp_int[p][CutIndex]
	 note number_conc "Number Concentration [p/cc]"
	 note number_conc "Integrated up to " + num2str(CutDiameter) + " nm"
	 
	 integrate/DIM=1/METH=1 dSdlogDp /X=logDp /D=dXdlogDp_int
	 surfacearea = dxdlogdp_int[p][CutIndex]
	 surfacearea = dxdlogdp_int[p][CutIndex]
	 note surfacearea "Surface area Concentration [nm^2/m3]"
	 note surfacearea "Integrated up to " + num2str(CutDiameter) + " nm"
	 
	 integrate/DIM=1/METH=1 dVdlogDp /X=logDp /D=dXdlogDp_int
	 volume_conc = dxdlogdp_int[p][CutIndex]
	 volume_conc = dxdlogdp_int[p][CutIndex]
	 note volume_conc "Volume Concentration [um^3/cm3]"
	 note volume_conc "Integrated up to " + num2str(CutDiameter) + " nm"
	 Mass_conc = volume_conc * density
	 note/K mass_conc "Assumes density of " + num2str(density) + " g/cm^3"
	 note mass_conc "In ug/m^3"
	 killwaves/z dXdlogDp_int, logDp
	 
End


//////////////////////////////////////////////////////////////////////////////////////////////////
// Function to average size distributions from image plots
// creats waves in the same folder as where the image plot information exist
Function AveDistFromImage(newname)
	string newname

	variable xposA = pcsr(A)
	variable xposB = pcsr(B)
	variable xpos_low, xpos_high
	xpos_low = min(xposA,xposB)
	xpos_high = max(xposA,xposB)
	variable npnts = xpos_high-xpos_low
	
	wave dXdlogDp = CsrWaveRef(A)
	string dXdlogDp_Str = CsrWave(A)
	variable V_smps = dimsize(dXdlogDp,1)
	
	string wvpath = GetWavesDataFolder(dxdlogdp,1)
	setdatafolder $wvpath	
	
	make/o/d/n=(V_smps) $newname=0, dXdlogDp_temp=0
	wave dXdlogDp_new = $newname
	make/o/d/n=(V_smps) $(newname+"_sdev")=0
	wave dXdlogDp_sdev = $(newname+"_sdev")
	
	note/k dXdlogDp_new "Averaged from " + dXdlogDp_Str + " wave, from " + num2str(xpos_low) + " to " + num2str(xpos_high)
	
	variable i, counter=0
	
	// average
	for(i=xpos_low;i<=xpos_high;i+=1)
		dxdlogdp_temp = dxdlogdp[i][p]
		wavestats/q dxdlogdp_temp
		if(V_numnans>4)
			// do nothing
		else
			dxdlogdp_new += dxdlogdp_temp
			counter+=1
		endif
	endfor	
	
	dxdlogdp_new/=counter
	
	// standard deviation
	for(i=xpos_low;i<=xpos_high;i+=1)
		dxdlogdp_temp = dxdlogdp[i][p]
		wavestats/q dxdlogdp_temp
		if(V_numnans>4)
			// do nothing
		else
			dxdlogdp_sdev += (dxdlogdp_temp-dxdlogdp_new)^2
			counter+=1
		endif
	endfor	
	
	dxdlogdp_sdev/=counter
	dxdlogdp_sdev = sqrt(dxdlogdp_sdev)
	
	killwaves/z dxdlogdp_temp
End

//**********************************************************************
Function RemoveBadFromSizeDist()
	
	string sizedist_str="dNdlogDp", twave_str="StartTime"
	prompt sizedist_str, "Size Distribution to Fix", popup WaveList("*",";","")
	prompt twave_str, "Time Wave", popup WaveList("*",";","")
	DoPrompt "Select Waves", twave_str, sizedist_str
	
	wave sizedist = $sizedist_str
	wave twave = $twave_str
	
	wave badstarttime = root:info:badperiodstartSMPS
	wave badstoptime = root:info:badperiodstopSMPS
	
	if(waveexists(badstarttime)==0 || waveexists(badstoptime)==0 || numpnts(twave) != dimsize(sizedist,0))
		abort "something is wrong."
	endif
	
	variable nb = numpnts(badstarttime)
	variable idex_start, idex_stop
	variable i
	
	for(i=0;i<nb;i+=1)
		findvalue/z /V=(badstarttime[i]) twave
		idex_start = V_value
		findvalue/z /V=(badstoptime[i]) twave
		idex_stop = V_value
		
		sizedist[idex_start,idex_stop][] = nan
	endfor
End

////*****************************************************************************************
//Function InterpSizeDistToNewDp(str_Dp_Old, str_dXdlogDp_Old, str_dp_new)
//	// uses linear interpolation
//	// sets values outside of size range of initial distribution to zero
//	string str_dp_old
//	string str_dXdlogDp_old
//	string str_dp_new
//
//	string str_dXdlogDp_new = str_dxdlogdp_old + "_interp"
//	
//	wave dp_old = $str_dp_old
//	wave dXdlogdp_old = $str_dxdlogdp_old
//	wave dp_new = $str_dp_new
//	
//	if(waveexists(dp_old)!=1 || waveexists(dxdlogdp_old)!=1 || waveexists(dp_new)!=1)
//		abort "one or more waves doesn't exist"
//	endif
//	
//	variable nd_old = numpnts(dp_old)
//	variable nd_new = numpnts(dp_new)
//	wavestats/q dp_old
//	variable old_max = V_max
//	variable old_min = V_min
//	
//	make/o/d/FREE/n=(nd_old) logDp_old, dX_old
//	make/o/d/FREE/n=(nd_new) logDp_new
//	make/o/d/n=(nd_new) $str_dxdlogdp_new
//	wave dxdlogdp_new = $str_dxdlogdp_new
//	
//	logDp_old = log(Dp_old)
//	logDp_new = log(Dp_new)
//	dX_old = dxdlogdp_old*logDp_old
//	
//	Interpolate2/T=1/N=200/I=3/F=0/Y=dx_new/X=dp_new dp_old, dx_old
//	wave dx_new
//	dxdlogdp_new = dx_new/logdp_new
//	dxdlogdp_new = dxdlogdp_new < 0 ? 0 : dxdlogdp_new
//	dxdlogdp_new = dp_new > old_max || dp_new < old_min ? 0 : dxdlogdp_new
//	
//	killwaves/z dx_new
////	return dxdlogdp_new
//End

//**********************************************************************************
function Bin_SizeDistByX(dXdlogDp,maskwavestr,binnedwavestr,binstart,binstop,nbins,[normalize])
	// bin size distribution data according to values in wave maskwave (e.g. photochemical age)
	// if there are any nans at a given time, make all sizes nan
	wave dXdlogDp	// size distribution wave
	string maskwavestr	// name of wave to use for masking
	string binnedwavestr	// name of new binned size distribution matrix
	variable binstart	// start of bin
	variable binstop	// last value of bin
	variable nbins	// number of bins
	variable normalize
	
	variable nRows=(dimsize(dXdlogDp,0))
	variable nDiam = dimsize(dXdlogDp,1)
	wave maskwave = $maskwavestr
	
	if(nRows != numpnts(maskwave))
		abort "Fool, your waves aren't the same length and they need to be!"
	endif
	if(waveexists(dXdlogDp)==0)
		Abort "You got yourself the wrong wave. Try again."
	endif
	if(ParamIsDefault(Normalize))
		Normalize = 0
	endif
	
	variable i, j
	variable binleft, binright
	variable binstep = (binstop-binstart)/nbins
	
	make/o/d/n=(nRows) dN
	make/o/d/n=(nBins,nDiam) $binnedwavestr
	make/o/d/n=(nRows,nBins) MatrixWaveTemp
	make/o/d/n=(nBins+1) BinX_im
	wave newwave = $binnedwavestr
	SetScale x, BinStart+BinStep/2, BinStop+BinStep/2, newwave//(BinStart+BinStep/2), BinStop+BinStep/2, newwave
	note/K newwave "Binned according to " + maskwavestr
	BinX_im = BinStart + binstep*x

	string str = binnedwavestr + "_sem"
	make/o/d/n=(nBins,nDiam) $str
	wave newwave_sem = $str
	string PointsPerBin = ""

	for(j=0;j<nDiam;j+=1)
		dN = dXdlogDp[p][j]	// get one column (i.e. size)
		for(i=0;i<nBins;i+=1)
			binleft = binstart+i*binstep
			binright = binstart+(i+1)*binstep
			extract dN, mywave, maskwave[x] >= binleft && maskwave[x] < binright // extract over bins
			if(numpnts(mywave)!=0)
				wavestats/q mywave
				newwave[i][j] = V_avg
				newwave_sem[i][j] = V_sem
			else
				newwave[i][j] = nan
				newwave_sem[i][j] = nan
				V_npnts = 0
			endif
			if(j==0)
				PointsPerBin += num2str(V_npnts)+";"
			endif
		endfor
	endfor
	note newwave "Points per bin is: " + PointsPerBin
	
	// Normalize size distributions to 1
	if(normalize==1)
		make/o/d/n=(nDiam)/FREE dNdlogDp_Single
		for(i=0;i<nBins;i+=1)
			dNdlogDp_single = newwave[i][p]
			wavestats/q dNdlogDp_single
			newwave[i][] = dNdlogDp_single[q]/V_max
		endfor
		note newwave "Each bin is normalized to max value"
	endif
	
	killwaves/z mywave, dN
end

//////////////////////////////////////////////////////////////////////////////////////////////////
// Function to determine the average size distribution between two points of size dist'n time series
// creates waves in the same folder as where the size dist'n information exist
// Assumes Rows = time, Columns = Diameter
// appends "_avg" to the wave
Function AveDistFromMatrix(wvstr,[start,stop,SetScaleStart, SetScaleDelta,PlusMinus,newwavestr,maskwavestr])
	string wvstr
	variable start
	variable stop
	variable SetScaleStart
	variable SetScaleDelta
	variable PlusMinus	// 0(default) = do nothing, 1 = create +/- 1 standard deviation
	string newwavestr
	string maskwavestr 
	
	wave dXdlogDp = $wvstr
	if(waveexists(dXdlogDp)==0)
		print "The wave upon which you wish to act does not exist. Thou must tryeth again."
	endif

	if(paramisdefault(newwavestr))
		newwavestr = wvstr + "_avg"
	endif
	
	variable nRows = dimsize(dXdlogDp,0)
	variable nDiam = dimsize(dXdlogDp,1)
	
	if(ParamIsDefault(start))
		start = 0
	endif
	if(ParamIsDefault(stop))
		stop = nRows
	endif
	if(ParamIsDefault(PlusMinus))
		PlusMinus=0
	endif
	if(ParamIsDefault(MaskWaveStr))
		MaskWaveStr = "noMask"
	else
		wave MaskWave = $MaskWaveStr
	endif

	variable npnts = stop-start
	
	if(PlusMinus==0)
		make/o/d/n=(nDiam) $newwavestr=0, dXdlogDp_temp=0
	else
		make/o/d/n=(nDiam,3) $newwavestr=0, dXdlogDp_temp=0
	endif
	wave dXdlogDp_new = $newwavestr
	
	variable i, counter=0
	
	if(stringmatch(maskwavestr,"noMask"))
		for(i=start;i<stop;i+=1)
			dxdlogdp_temp = dxdlogdp[i][p]
			wavestats/q dxdlogdp_temp
			if(V_numnans>0)
				// do nothing
			else
				dxdlogdp_new += dxdlogdp_temp
				counter+=1
			endif
		endfor
		dxdlogdp_new/=counter
	else
		make/o/d/n=(nRows)/FREE OneSize
		for(i=0;i<nDiam;i+=1)
			OneSize = dxdlogdp[p][i]
			extract/FREE/O OneSize, Extracted, MaskWave==1
			wavestats/q Extracted
			dXdlogDp_new[i] = V_avg
		endfor
	endif
	
	
	if(PlusMinus==1)
		make/o/d/n=(nRows)/FREE OneSize
		for(i=0;i<nDiam;i+=1)
			OneSize = dxdlogdp[p][i]
			wavestats/q/r=(start,stop) OneSize
			dxdlogdp_new[i][0] = V_avg
			dxdlogdp_new[i][1] = V_avg-V_sem
			dxdlogdp_new[i][2] = V_avg+V_sem
		endfor
	endif
	
	if(ParamIsDefault(SetScaleStart)==0 && ParamIsDefault(SetScaleDelta)==0)
		setscale/P x, SetScaleStart, SetScaleDelta, dxdlogdp_new
	endif
	
	killwaves/z dxdlogdp_temp
End

//**********************************************************************************
Function NormalizeEachRowOfDistn(dXdlogDp_str,[Pts2Avg])
	// normalize a 2D size distribution matrix (rows = time, columns = diameter)
	// to max value. Well, actually to average of Pts2Avg highest points.
	string dxdlogdp_str
	variable Pts2avg
	
	wave dxdlogdp = $dxdlogdp_str
	if(waveexists(dXdlogDp)==0)
		abort "The wave upon which you wish to act does not exist. Thou must tryeth again."
	endif
	
	if(ParamIsDefault(Pts2avg))
		pts2avg = 3
	endif
	
	string normstr = dxdlogdp_str + "_norm"
	variable nRows = dimsize(dxdlogdp,0)
	variable nDiam = dimsize(dxdlogdp,1)
	make/o/d/n=(nDiam)/FREE tempwave, sortedwave
	make/o/d/n=(nRows,nDiam) $normstr
	wave normwave = $normstr
	note normwave "Normalized data associated with " + dxdlogdp_str
	note normwave "Normalized to average of max " + num2str(3) + " points"
	variable i
	
	for(i=0;i<nRows;i+=1)
		tempwave = dxdlogdp[i][p]
		sortedwave = tempwave
		sort sortedwave, sortedwave
		wavestats/q/R=[nDiam-pts2avg,nDiam-1] sortedwave
		normwave[i][] = tempwave[q]/V_avg
	endfor	
End

//*******************************************************************************************
Function CleanUpSizeDist_tseries([dNdlogDpStr,dpStr,SMPSorAPS,LL,UL])
	// This function can be used to remove "bad" periods of data, for example when a filter was on
	// or when there was some crazy spike in the signal that is not real
	string dNdlogDpStr
	string dpStr
	string SMPSorAPS	// "SMPS" or "APS"
	variable LL	// lower limit integrated number concentration
	variable UL	// upper limit single point
	
	if(ParamIsDefault(SMPSorAPS))
		SMPSorAPS = "SMPS"
	endif
	if(ParamIsDefault(LL))
		LL = 10
	endif
	if(ParamIsDefault(UL))
		UL = 1e5
	endif
	
	string dfMerged = "root:"+SMPSorAPS+":merged"
	setdatafolder $dfMerged
	
	if(ParamIsDefault(dNdlogDpStr))
		dNdlogDpStr = "dNdlogDp_"+SMPSorAPS
	endif
	if(ParamIsDefault(dpStr))
		dpStr = "dia_SMPS"
	endif
	
	wave dNdlogDp = $dNdlogDpStr
	wave dp = $dpStr
	variable nRows = dimsize(dNdlogDp,0)
	variable nCols = dimsize(dNdlogDp,1)
	make/o/d/n=(nCols)/FREE logDp
	logDp = log(dp)
	make/o/d/n=(nRows) Concentration
	
	integrate/T/DIM=1 dNdlogDp /X=logDp /D=IntWave
	Concentration = IntWave[p][nCols-1]
	
	dNdlogDp = Concentration[p] < LL || Concentration[p] > UL ? nan : dNdlogDp[p][q]
	
	integrate/T/DIM=1 dNdlogDp /X=logDp /D=IntWave
	Concentration = IntWave[p][nCols-1]
	KillWaves/Z IntWave, Concentration
End

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To fit a time-series of size distributions to a bimodal log-normal
function fit_2modeLogNorm_tseries(dXdlogDpStr,DiameterStr,low_size, high_size)
	string dXdlogDpStr
	string DiameterStr
	variable low_size, high_size
		
	wave dXdlogDp = $(dXdlogDpStr)
	wave dia = $(DiameterStr)
	variable npnts_Time = dimsize(dxdlogdp,0)
	variable npnts_Diam = dimsize(dXdlogDp,1)
		
	make/o/n=(npnts_Time,6) $("dXdlogDp_fit_coefs")
	wave int_v_fit = $("smps_int_V_FIT")
	Make/O/T/N=7 w_constraints
	wave/t w_constraints = w_constraints
	w_constraints[0] = {"K2 > 20","K2 < 300","K3 > 0.3","K3 < 1.5","K5 > 100","K5 < 400","K6 > 0.3","K6 < 1.5"}	//"K0 >(-0.1)","K0 < 10",
	make/o/n=(npnts_Time,7) $("fit_coefs")
	wave fit_coefs = $("fit_coefs")
	variable i,j
	variable scale

	duplicate/o dia, $("logDp")
	wave log_dia = logDp
	log_dia = log(dia)
	
	make/o/n = (npnts_Time) fit_int_V
	wave fit_int_V
	make/o/n=1000 fit_dia, fit_log_dia
	make/o/n=(npnts_Time,1000) Fit_Distr
	wave fit_dia, fit_log_dia
	fit_dia = low_size+(high_size-low_size)/1000*p
	setscale/P y, low_size, (high_size-low_size)/1000, fit_distr
	fit_log_dia = log(fit_dia)
	
	for (i=0;i<npnts_Time;i+=1)
		if (numtype(dXdlogDp[i][0])==0)
			make/o/n=(npnts_Diam) SingleDist
			wave SingleDist
			SingleDist = dXdlogDp[i][p]
			wavestats/q SingleDist
			scale = V_max
			make/o/n=7 w_coefs = {0,V_max/1.5,Dia[V_maxloc]-20,0.5,V_max/1.5,Dia[V_maxloc]+20,0.5}
			wave w_coefs
			FuncFit/Q/H="1000000"/N LogNormal_2modes kwCWave=w_coefs,  SingleDist /X=dia /D /C=w_constraints ///D=newfit 
			wave/z fit_SingleDist

			fit_distr[i][] = w_coefs[0]+w_coefs[1]*exp(-(ln(fit_dia[q]/w_coefs[2])/w_coefs[3])^2)+w_coefs[4]*exp(-(ln(fit_dia[q]/w_coefs[5])/w_coefs[6])^2)
			fit_coefs[i][] = w_coefs[q]
		endif
	endfor
	
	killwaves/Z w_coef, w_constraints, w_paramconfidenceinterval, w_sigma
	killwaves/Z SingleDist, fit_SingleDist, fit_source2
	Killwaves/Z logDp, smps_V_fit_coefs
	killwaves/z w_coefs, fit_int_V
end

//Function InterpSizeDist_toNewDp(dNdlogDp,DpOld,DpNew,NewName)
//	// Function to convert  one size distribution (or a time-series of size distributions) from one size range to another
//	// works through integration, interpolation and differentiation
//	// the differentiation step seems to impart some smoothing to the data, quite annoyingly
//	// written by CD Cappa on 12/31/13
//	wave dNdlogDp	// dN/dlogDp for original size distribution
//	wave DpOld	// diameters of original size distribution
//	wave DpNew // diameters of new size distribution
//	string NewName // name of new size distribution
//	
//	variable NumCols = dimsize(dNdlogDp,1)
//	variable NumRows = dimsize(dNdlogDp,0)
//	variable numDpOld = numpnts(DpOld)
//	variable numDpNew = numpnts(DpNew)
//	variable multidimensional
//	variable numDp
//	variable i
//	if(NumCols>1)
//		multidimensional = 1
//		numDp = NumCols
//	else
//		multidimensional = 0
//		numDp = NumRows
//	endif
//	if(multidimensional==1 && numDpOld != NumCols)
//		abort "dNdlogDp is multidimensional and number of columns doesn't match DpOld"
//	elseif(multidimensional==0 && numDpOld != NumRows)
//		abort "dNdlogDp is a single size distribution and number of rows doesn't match DpOld"
//	endif
//	
//	wavestats/q DpOld
//	variable minDp_Old = V_min
//	variable maxDp_Old = V_max
//	
//	wavestats/q DpNew
//	variable minDp_new = V_min
//	variable maxDp_new = V_max
//	
//	make/o/d/n=(numDpOld) logDp_Old = log(DpOld)
//	make/o/d/n=(numDpNew) logDp_New = log(DpNew)
//		
//	if(wavetype(dNdlogDp)==2)
//		redimension/S logDp_old
//		redimension/S logDp_new
//	endif
//	
//	if(multidimensional==0)
//		Integrate/METH=1 dNdlogDp /X=logDp_Old /D=dNdlogDp_old_INT
//	else
//		Integrate/DIM=1/METH=1 dNdlogDp /X=logDp_Old /D=dNdlogDp_old_INT
//	endif
//	
//	if(multidimensional==0)
//		make/o/d/n=(numDpNew) dNdlogDp_New
//		make/o/d/n=(numDpNew) dNdlogDp_New_INT
//	else
//		make/o/d/n=(numRows,numDpNew) dNdlogDp_New = nan
//		make/o/d/n=(numRows,numDpNew) dNdlogDp_New_INT = nan
//		make/o/d/n=(numDpNew) dNdlogDp_new_single_INT
//		make/o/d/n=(numDpOld) dNdlogDp_old_single_INT
//	endif
//	
//	if(multidimensional==0)
//		wavestats/q dNdlogDp_old_int
//		if(V_numnans==0)
//			Interpolate2/T=1/I=3/Y=dNdlogDp_new_INT/X=logDp_New logDp_Old, dNdlogDp_old_INT
//			Differentiate dNdlogDp_new_INT/X=logDp_New/D=dNdlogDp_new
//		else
//			dNdlogDp_new = nan
//		endif
//	else
//		for(i=0;i<numRows;i+=1)
//			dNdlogDp_old_single_int = dNdlogDp_Old_INT[i][p]
//			wavestats/q dNdlogDp_old_single_int
//			if(V_numnans==0)
//				Interpolate2/T=2/I=3/Y=dNdlogDp_new_single_INT/X=logDp_New logDp_Old, dNdlogDp_old_single_INT
//				dNdlogDp_new_INT[i][] = dNdlogDp_new_single_INT[q]
//				Differentiate dNdlogDp_new_single_INT/X=logDp_New/D=dNdlogDp_new_single_INT_DIF
//				dNdlogDp_new_single_INT_DIF = DpNew < minDp_Old || DpNew > maxDp_Old ? nan : dNdlogDp_new_single_INT_DIF
//				dNdlogDp_New[i][] = dNdlogDp_new_single_INT_DIF[q]
//			else
//				dNdlogDp_New[i][] = nan
//			endif
//		endfor
//		
//	endif
//	
//	if(multidimensional==0 && stringmatch(NewName,"dNdlogDp_New")==0)
//		make/o/d/n=(numDpNew) $NewName
//		wave NewWave = $NewName
//		NewWave = nan
//		NewWave = dNdlogDp_new
//		killwaves/z dNdlogDp_new
//	elseif(multidimensional==1 && stringmatch(NewName,"dNdlogDp_New")==0)
//		make/o/d/n=(numRows,numDpNew) $NewName
//		wave NewWave = $NewName
//		NewWave = nan
//		NewWave = dNdlogDp_new
//		killwaves/z dNdlogDp_new
//	endif
//	
//	killwaves/z dNdlogDp_new_single_INT, dNdlogDp_old_single_INT
//	killwaves/z DpNew_expanded, dNdlogDp_new_INT, dNdlogDp_new_Single_INT
//	killwaves/z dNdlogDp_old_single_INT, dNdlogDp_Old_INT
//	killwaves/z dNdlogDp_new_single_INT_DIF
//	killwaves/z logDp_Old, logDp_New
//	killwaves/z logDp_new, dNdlogDp_old_INT
//
//End

//************************************************************************************
Function Interp_dNdlogDp_toNewDp(dNdlogDp,DpOld,DpNew,NewName)
	// Function to convert  one size distribution from one size range to another
	// may run into problems if new size range is outside original size range
	// Written by CDC on 9/28/14
	wave dNdlogDp	// dN/dlogDp for original size distribution
	wave DpOld	// diameters of original size distribution
	wave DpNew // diameters of new size distribution
	string NewName // name of new size distribution
	
	variable nbins_old = numpnts(dpOld)
	variable nbins_new = numpnts(DpNew)
	variable i
	
	make/o/d/n=(nbins_new)/FREE Dp_Low, Dp_High, Index_Low, Index_High
	make/o/d/n=(nbins_new)/FREE logDpNew = log(DpNew)
	make/o/d/n=(nbins_old)/FREE logDpOld = log(DpOld)
	make/o/d/n=(nbins_new) $NewName
	wave dNdlogDp_New = $NewName
	
	Differentiate/METH=0 logDpNew /D=dlogDpNew
	Dp_Low = 10^(logDpNew - 0.5*dlogDpNew)
	Dp_high = 10^(logDpNew + 0.5*dlogDpNew)
	
	for(i=0;i<nbins_new;i+=1)
		FindLevel/Q DpOld, Dp_Low[i]
		If(V_flag==0)
			Index_Low[i] = V_levelX
		else
			Index_Low[i] = 0
		endif
		FindLevel/Q DpOld, Dp_High[i]
		if(V_flag==0)
			Index_high[i] = V_levelX
		else
			Index_high[i] = nbins_old-1
		endif
	endfor
	
	Integrate/METH=1 dNdlogDp /X=logDpOld /D=dN
	for(i=0;i<nbins_new;i+=1)
		dNdlogDp_new[i] = (dN[Index_high[i]] - dN[Index_low[i]])/dlogDpNew[i]
	endfor
	
	KillWaves/z dN, dlogDpNew
end


//*********************************************************************************************
Function Interp_dNdlogDp_toNewDp_2D(dNdlogDp,DpOld,DpNew,NewName)
	// Function to convert  one size distribution from one size range to another
	// may run into problems if new size range is outside original size range
	// Written by CDC on 9/28/14
	// updated 05/09/16 to use differentiation for reverse step
	wave dNdlogDp	// dN/dlogDp for original size distribution, 2D matrix with time as rows, diameter as columns
	wave DpOld	// diameters of original size distribution
	wave DpNew // diameters of new size distribution
	string NewName // name of new size distribution
	
	variable nrows = dimsize(dNdlogDp,0)
	variable nbins_old = numpnts(dpOld)
	variable nbins_new = numpnts(DpNew)
	variable i
	
	make/o/d/n=(nbins_new)/FREE Dp_Low, Dp_High, Index_Low, Index_High
	make/o/d/n=(nbins_new) logDpNew = log(DpNew)
	make/o/d/n=(nbins_old)/FREE logDpOld = log(DpOld)
	make/o/d/n=(nrows,nbins_new) $NewName
	wave dNdlogDp_New = $NewName
	dNdlogDp_new = nan
	make/o/d/n=(nrows,nbins_new) dN_new
	
	Differentiate/METH=0 logDpOld /D=dlogDpOld
	Differentiate/METH=0 logDpNew /D=dlogDpNew
	Integrate/DIM=1/METH=1 dNdlogDp /X=logDpOld /D=dN
	make/o/d/n=(nbins_old) dN_1row
	for(i=0;i<nrows;i+=1)
		dN_1row = dN[i][p]
		wavestats/q dN_1row
		if(V_numnans > 30)
			// do nothing
			dN_new[i][] = nan
		else
			Interpolate2/T=2/E=2/I=3/Y=dN_1row_interp/X=dpnew dpold, dN_1row
			dN_new[i][] = dN_1row_interp[q]
		endif	
	endfor
	
	Differentiate/DIM=1 dN_new/D=dN_new_DIF // added 05/09/16
	dNdlogDp_new = dN_new_DIF/dlogDpNew[q]
	
//	for(i=nbins_new;i>0;i-=1) // removed 05/09/16
//		dNdlogDp_new[][i] = (abs(dN_new[p][i]-dN_new[p][i-1]))/dlogDpNew[i]
//	endfor
	
	dNdlogDp_new[][] = dpNew[q] < dpold[0] ? nan : dNdlogDp_new[p][q]
	dNdlogDp_new[][] = dpNew[q] > dpold[nbins_old-1] ? nan : dNdlogDp_new[p][q]
	for(i=0;i<nbins_new;i+=1)
		dNdlogDp_new[i][] = numtype(dNdlogDp_new[i][q])==2 ? 0 : dNdlogDp_new[i][q] // updated 05/09/16
	endfor

	KillWaves/z dN, dlogDpNew, dN_1row, dN_1row_interp, dN_new, dlogDpOld, logDpNew
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function DisplaySizeDist_tseries(wv,dp,start,stop,normalize)
	wave wv
	wave dp
	variable start,stop, normalize
	
	display
	variable i
	string norm_str
	for(i=start;i<=stop;i+=1)
		appendtograph wv[i][] vs dp
	endfor
	ModifyGraph log(bottom)=1
	ModifyGraph lsize=2
	SetAxis left 0,1
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function NormalizeSizeDist_tseries(wvstr,[LL,UL])
	string wvstr
	variable UL, LL
	
	wave wv = $wvstr
	variable npnts = dimsize(wv,1)
	variable ntime = dimsize(wv,0)
	if(paramisdefault(LL))
		LL = 0
	endif
	if(paramisdefault(UL))
		UL = npnts
	endif
	
	variable i
	make/o/d/n=(npnts)/FREE single
	make/o/d/n=(ntime,npnts) $(wvstr+"_norm") = nan
	wave wvnorm = $(wvstr+"_norm")
	for(i=0;i<ntime;i+=1)
		single = wv[i][p]
		wavestats/q/R=[LL,UL] single
		wvnorm[i][] = single[q]/V_max
	endfor
end

///////////////////////////////////////// SEMS BELOW ///////////////////////////////////////////////////////////////
Function Load_SEMS_files()
	// CDC 07/05/17: changed so that data are not loaded in root folder

	setdatafolder root:
	
	string extension = "dat"
	Variable refNum
	String message = "Select one or more CONC files"
	String/G file_list
	String fileFilters = "Data Files (*."+extension+"):."+extension+";"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	file_list = S_fileName
	if (strlen(file_list) == 0)
		print "==================================="
		print "No Data Files Chosen for Loading.  Aborting."
		print "==================================="
		Close/A
		setdatafolder root:
		abort
	endif
	Close/A
	file_list = replacestring("\r",file_list,";")
	file_list = sortlist(file_list,";",16)
	// items in list will have names like "SEMS_CONC_XXXXXX_XXXXXX"

	variable n_files = itemsinlist(file_list,";")
	variable i, j
	
	string SEMS_base_df = "root:SEMS"
	string SEMS_loadDF = "root:SEMS:loading"
	string SEMS_df
	string df_str
	variable df_items
	if(datafolderexists(SEMS_base_df)==0)
		newdatafolder $SEMS_base_df
	endif
	
	if(datafolderexists(SEMS_loadDF)==0)
		newdatafolder $SEMS_loadDF
	endif	

	for(i=0;i<n_files;i+=1) // loop over files
// CDC 07/05/17		setdatafolder root:
		setdatafolder $SEMS_loadDF // CDC 07/05/17
		killwaves/a/z
		
		semsloader(file_list,i,loadDF=SEMS_loadDF) // this is a brechtel procedure
		
		df_str = stringfromlist(i,file_list,";")
		df_items = itemsinlist(df_str,":")
		df_str = stringfromlist(df_items-1,df_str,":")
		df_str = removefromlist("SEMS_CONC",df_str,"_")
		df_str = stringfromlist(0,df_str,".")
		df_str = "x" + df_str

		setdatafolder $SEMS_base_df
		SEMS_df = SEMS_base_df + ":" + df_str
		if(datafolderexists(SEMS_df)==0)
			newdatafolder $SEMS_df
		endif
// CDC		setdatafolder root:
		setdatafolder $SEMS_loadDF
		
		wave starttime
		wave startdate
		wave startsec, endsec
		variable n_time = numpnts(starttime)
		variable yr, da, mo
		string datestr
	
		// create a useful timewave
		make/o/d/n=(n_time) $(SEMS_df+":timestart")
		wave timestart =$(SEMS_df+":timestart") 
		setscale/p d, 0, 1, "dat", timestart
		// time conversion updated 6/24/17 to deal with data that goes over midnight
		datestr = num2istr(startdate[0])
		yr = 2000 + str2num(datestr[0,1])
		mo = str2num(datestr[2,3])
		da = str2num(datestr[4,5])
		if(startsec[0] > endsec[0])
			timestart = date2secs(yr,mo,da) + starttime[0] + endsec - endsec[0]
		else
			timestart = date2secs(yr,mo,da) + starttime[0] + startsec - startsec[0]
		endif

		// get image plot time
		make/o/d/n=(n_time+1) $(SEMS_df + ":timestart_im") = 0
		wave timestart_im = $(SEMS_df + ":timestart_im")
		setscale/p d, 0, 1, "dat", timestart_im
		make/o/d/n=(n_time+1) dt
		dt[0,n_time-2] = timestart[x+1]-timestart[x]
		dt[n_time-1,] = dt[n_time-2]
		timestart_im[0] = timestart[0] - dt[0]/2
		for(j=0;j<n_time;j+=1)
			timestart_im[j+1] = timestart_im[j] + dt[j]
		endfor
		
		// Get average diameter across all scans
		// Why? because fuck using an individual diameter wave for every single scan
		wave nm_scan1
		variable n_diam = numpnts(nm_scan1)
		make/o/d/n=(n_diam) $(SEMS_df + ":Dpm") = 0
		wave Dpm = $(SEMS_df + ":Dpm")
		
		for(j=0;j<n_time;j+=1)
			wave dp_single = $("nm_scan"+num2istr(j+1))
			Dpm += dp_single
		endfor
		Dpm /= n_time
		// get image plot diameter
		make/o/d/n=(n_diam+1) $(SEMS_df + ":Dpm_im") = 0
		wave Dpm_im = $(SEMS_df + ":Dpm_im")
		make/o/d/n=(n_diam+1) logDp_temp, dlogDp_temp
		logDp_temp[0,n_diam-1] = log(Dpm)
		logDp_temp[n_diam,] = logDp_temp[n_diam-1]
		dlogdp_temp[0,n_diam-1] = logdp_temp[x+1]-logDp_temp[x]
		dlogdp_temp[n_diam-1,] = dlogdp_temp[n_diam-2]
		Dpm_im[0] = 10^(log(Dpm[0]) - dlogDp_temp[0]/2)
		for(j=0;j<n_diam;j+=1)
			Dpm_im[j+1] = 10^(log(Dpm_im[j]) + dlogdp_temp[j])
		endfor
			
		//Get dNdlogDp
		make/o/d/n=(n_time,n_diam) $(SEMS_df + ":dNdlogDp") = 0
		wave dNdlogDp = $(SEMS_df + ":dNdlogDp")
		
		for(j=0;j<n_time;j+=1)
			wave dndlogdp_single = $("dndlogdp_scan"+num2istr(j+1))
			dndlogdp[j][] = dndlogdp_single[q]
		endfor
		
		// Get Sample Flow
		make/o/d/n=(n_time) $(SEMS_df+":SampleFlow")
		wave sampleflow =$(SEMS_df+":SampleFlow")
		wave SF = $(SEMS_loadDF+":sampavg") // name of the sample flow wave
		SampleFlow = SF
		note/K SampleFlow "Sample flow rate measured by the LFE; units = lpm"

	endfor
// CDC	setdatafolder root:
	setdatafolder $SEMS_loadDF
	killwaves/a/z
	setdatafolder $SEMS_base_DF
	killdatafolder/z $SEMS_loadDF
	
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to merge multiple SEMS files into one long time series
// Currently not general b/c must set folders to be used (must update "fldrs_str" wave)
Function MergeSEMS_tseries(datatype,deleteOldFolders)
	string datatype // e.g. dNdlogDp, dN, dNdDp, etc. 
	variable deleteOldFolders // 0 = no, 1 = yes

	string df_SEMS = "root:SEMS:"
	setdatafolder $df_SEMS
	string df_SEMS_Merged = "root:SEMS:Merged"
	if(datafolderexists(df_SEMS_Merged)==0)
		newdatafolder $df_SEMS_Merged
	endif
		
	variable ndf = countobjects(":",4)	//directories with files
	string/g datafolder_list,datafolder_str, file_name_str
	//get datafolder list, check for "x" at begining due to start of file with number, sort alpha numerically
	datafolder_list = stringbykey("FOLDERS",DataFolderDir(1))
	//remove ("Merged") datafolder name from list
	if (stringmatch(datafolder_list,"*Merged*"))
		datafolder_list = RemoveFromList("Merged",datafolder_list,",")
		ndf -= 1
	endif
	//sort datafolders alpha numerically.  NOTE this will not sort mixed "x" directories with other directories well...
	datafolder_list = sortlist(datafolder_list,",",16)	

	setdatafolder df_SEMS_merged

	make/o/n=(ndf) wv_npnts = nan
	
	variable V_npntstot 
	variable V_dia
	
	string cdf_str
	string wvstr_dia1 = "dia_"
	string wvstr_time = "dt_"
	
	string str_dndlogdp = "dNdlogDp"
	string str_time = "timestart"
	string str_time_stop = "timestop"
	string str_dia = "Dpm"
	string str_SF = "SampleFlow" // CDC 07/05/17

	string current_str
	
	variable i, j
		
	make/o/d/n=(V_dia) Dpm = nan //dia_SEMS = nan
	wave dia_SEMS = Dpm
	
	make/o/d/n=(V_dia+1) Dpm_im = nan //dia_SEMS_im = nan
	wave dia_SEMS_im = Dpm_im
	
	string dt_list = "", dt_im_list="", sdm_list="", SF_list=""
	
	for(i=0;i<ndf;i+=1)
		cdf_str = df_SEMS + stringfromlist(i,datafolder_list,",")//[i]
		setdatafolder $cdf_str
		
		current_str = "Dpm" // wvstr_dia1 + stringfromlist(i,datafolder_list,",")
		wave current = $current_str
		if(i==0) 
			duplicate/o current, $(df_SEMS+"merged:Dpm")//dia_SEMS")//dia_SEMS = current
		endif
		
		current_str = str_dia + "_im" // "dia_" + stringfromlist(i,datafolder_list,",") + "_im"
		wave current = $current_str
		if(i==0)
			duplicate/o current, $(df_SEMS+"merged:Dpm_im")//dia_SEMS_im")//dia_SEMS_im = current
		endif
		
		dt_list += df_SEMS+stringfromlist(i,datafolder_list,",")+":"
		dt_list += str_time + ";" //wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		
		if(i<ndf-1)
			dt_im_list += df_SEMS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += str_time + ";" // wvstr_time+stringfromlist(i,datafolder_list,",")+";"
		else
			dt_im_list += df_SEMS+stringfromlist(i,datafolder_list,",")+":"
			dt_im_list += str_time + "_im;" //wvstr_time+stringfromlist(i,datafolder_list,",")+"_im;"
		endif
		
		sdm_list += df_SEMS+stringfromlist(i,datafolder_list,",")+":"
		sdm_list += str_dndlogdp + ";" //"sdm_"+stringfromlist(i,datafolder_list,",")+";"
		
		SF_list += df_SEMS+stringfromlist(i,datafolder_list,",")+":"
		SF_list += str_SF + ";" //wvstr_time+stringfromlist(i,datafolder_list,",")+";"

	endfor
	
	setdatafolder df_SEMS_Merged	
	concatenate/o dt_list, dt_SEMS
	redimension/n=(-1,0) dt_SEMS
	setscale/p d 0,1, "dat", dt_SEMS
	duplicate/o dt_SEMS $str_time
	killwaves/z dt_SEMS	
	
	concatenate/o dt_im_list, dt_SEMS_im
	redimension/n=(-1,0) dt_SEMS_im
	setscale/p d 0,1, "dat", dt_SEMS_im
	duplicate/o dt_SEMS_im $(str_time+"_im")
	killwaves/z dt_SEMS_im
	
	concatenate/o/np=0 sdm_list, sdm_SEMS
	duplicate/o sdm_SEMS, $str_dndlogdp // $(datatype+"_SEMS")
	wave sdm_merged = $str_dndlogdp // $(datatype+"_SEMS")
	killwaves sdm_SEMS
	
	sdm_merged = numtype(sdm_merged)==2 ? 0 : sdm_merged
	
	concatenate/o SF_list, SF_SEMS
	redimension/n=(-1,0) SF_SEMS
	duplicate/o SF_SEMS $str_SF
	killwaves/z SF_SEMS
	
	killwaves/z wv_npnts
	
	if(deleteOldFolders==1)
		for(i=0;i<ndf;i+=1)
			cdf_str = df_SEMS + stringfromlist(i,datafolder_list,",")//[i]
			setdatafolder $df_SEMS
			killdatafolder/z $cdf_str
		endfor
	endif
	
	setdatafolder df_SEMS_Merged	
End


//***************************************************************************
Function SEMS_CorrectForFlow([SF_Reference])
	// run this only after running merging and run only once
	variable SF_reference
	
	if(paramisdefault(SF_reference))
		SF_reference = 0.37 // lpm
	endif

	string df_merged = "root:SEMS:Merged"
	string str_SF = "SampleFlow"
	string str_dN = "dNdlogDp"
	
	setdatafolder $df_merged
	wave SF = $str_SF
	wave dn = $str_dN
	
	make/o/d/n=(numpnts(SF)) SF_Correction
	SF_Correction = SF_Reference/SF
	smooth 4, SF_Correction
	
	dn[][] = dn[p][q] * SF_Correction[p]
	note/k dn "Corrected for flow variations using reference flow = " + num2str(SF_reference) + " lpm"

End

//////////////////////////////////////////////////////////////////////////////////////////////////
Function SEMSloader(file_list,filenum,[loadDF])
	string file_list
	variable filenum
	string loadDF // folder where to load things

//	setdatafolder root:
	setdatafolder $loadDF
	killwaves/a/z
		
	variable i					//loop counter
	variable j					//loop counter
	
//Filename Vars
	String invdata = ""
	invdata = stringfromlist(filenum,file_list,";")
	String rawdata = ""
	rawdata = replacestring("CONC",invdata,"RAW")
	Variable refnum

	String pathstr = ""
	variable pathitems
	pathitems = itemsinlist(invdata,":")
	pathstr = stringfromlist(pathitems-1,invdata,":")
	pathstr = removefromlist(pathstr,invdata,":")
	NewPath/C/Q/O DataPath pathstr
		
	String buffer = "" 			//line read buffer
	
	variable lineNumber, len 	//file length counters
	variable temp0
	variable temp1
	
//Scan Data
	variable numscans		//num scans from inv data file
	variable numscansraw		//num scans from raw data file
	variable numwaves		//waves loaded from inv data file
	variable numwavesraw		//waves loaded from raw data file
	variable numbincols		//number of bins from inv data file
	variable numbincolsraw	//number of bins from raw data file

//Temp wave name vars
	string s_info=""	
	string wavestr =""
	string wavestr1 =""
	string binstr=""
	string binstr1=""
	string binstr2=""
	string scanstr=""
	string scanstr1=""
	variable wavenum
	variable wavenum1
	
//-----------------------------------
//Read Inverted Data
//-----------------------------------
//	print "Select SEMS_CONC...dat file"
	
	string extension = "dat"
	String message = "Select one or more files"
	String fileFilters = "Data Files (*."+extension+"):."+extension+";"
	
//	Open/R /F=filefilters refNum as ""				// Display dialog
	Open/R /F=filefilters /P=DataPath refnum as invdata
	invdata = S_filename
	rawdata = replacestring("CONC",invdata,"RAW")
	pathitems = itemsinlist(rawdata,":")
	pathstr = stringfromlist(pathitems-1,invdata,":")
	rawdata = stringfromlist(pathitems-1,rawdata,":")
	pathstr = removefromlist(pathstr,invdata,":")
	NewPath/C/Q/O DataPath pathstr
	
//	print invdata
	if (refNum == 0)
		return -1						// User canceled
	endif

	lineNumber = 0
	
	do
			FReadLine refNum, buffer
			len = strlen(buffer)
			if (len == 0)
				break						// No more lines to be read
			endif
		lineNumber += 1
		
	while (1)
	
	Close refnum
	
	numscans = ( linenumber - 27 )
	//print "Num Scans"
	//print numscans
	
	//----------------------------------------------------
	
	loadwave /q/a/J/L={0,27,0,0,0}/k=0 invdata
	numwaves = V_flag
	numbincols = numwaves - 10
	
	rename wave0, startdate
	wave startdate
	rename wave1, starttime
	wave starttime
	rename wave2, endate
	wave endate
	rename wave3, endtime
	wave endtime
	rename wave4, startsec
	wave startsec
	rename wave5, endsec
	wave endsec
	rename wave6, scantemp
	wave scantemp
	rename wave7, scanpress
	wave scanpress
	rename wave8, scanbins
	wave scanbins
		
	for( i = 0; i < scanbins[0] + 1; i += 1 )		//rename bin conc autonamed waves
		wavenum = 9 + i
		wavestr = "wave" + num2str(wavenum)
		wave w = $wavestr
		
		if( i < scanbins[0] )
			binstr = "bin" + num2str( i + 1 )
			rename w, $binstr
		
		elseif( i == scanbins[0] )
			rename w, numflag
			wave numflag
		
		endif		
	endfor
	//edit
	for( i = 0; i < numscans  ; i += 1 )			//reorg binconc[scan] waves into scan[bin] waves
		scanstr = "dndlogdp_scan" + num2str( i + 1)
		make /o/n =(scanbins[i]) $scanstr
		wave scan = $scanstr
		//appendtotable $scanstr
		for( j = 0; j < scanbins[i]; j += 1 )
			binstr = "bin" + num2str( j + 1 )
			wave bin = $binstr
			scan[j] = bin[i]
			
		endfor	
	endfor
	
	
//-------------------------
//Read Raw File
//-------------------------
//	print "Select SEMS_RAW...dat file"
//	Open/R /F=filefilters refNum as ""				// Display dialog
	Open/R /F=filefilters /P=DataPath refNum as rawdata				// Display dialog
	rawdata = S_filename
//	print rawdata
	if (refNum == 0)
		return -1						// User canceled
	endif

	lineNumber = 0
	
	do
			FReadLine refNum, buffer
			len = strlen(buffer)
			if (len == 0)
				break						// No more lines to be read
			endif
		lineNumber += 1
		
	while (1)
	
	Close refnum
	
	numscansraw = ( linenumber - 22 )
	
	//print "Num Scans"
	//print numscansraw
	
	//----------------------------------------------------
	
	loadwave /q/a/J/L={0,21,0,0,0}/k=0 rawdata
	numwavesraw = V_flag
	numbincolsraw = numwavesraw - 10
	
	rename wave0, startdateraw
	wave startdateraw
	rename wave1, starttimeraw
	wave starttimeraw
	rename wave2, endateraw
	wave endateraw
	rename wave3, endtimeraw
	wave endtimeraw
	rename wave4, startsecraw
	wave startsecraw
	rename wave5, endsecraw
	wave endsecraw
	rename wave6, scantempraw
	wave scantempraw
	rename wave7, scanpressraw
	wave scanpressraw
	rename wave8, sampavg//raw
	wave sampavg//raw
	rename wave9, shavg
	wave shavg
	rename wave10, samprh
	wave samprh
	rename wave11, shrh
	wave shrh
	rename wave12, scanbinsraw
	wave scanbinsraw
		
//	for( i = 0; i < scanbinsraw[0] + 8; i += 1 )		//rename bin count autonamed waves
	for( i = 0; i < scanbinsraw[0]; i += 1 )		//rename bin count autonamed waves

		wavenum = 13 + i
		wavestr = "wave" + num2str(wavenum)
		wave w0 = $wavestr
		wavestr1 = "wave" + num2str(wavenum + scanbinsraw[0])
		wave w = $wavestr1
		
		if( i < scanbinsraw[0] )
			binstr = "bindia" + num2str( i + 1 )
			rename w0, $binstr
			binstr1 = "bincnt" + num2str( i + 1 )
			rename w, $binstr1
			
		elseif( i == scanbinsraw[0] + 0 )
			rename w, sattemp
			wave sattemp
			
		elseif( i == scanbinsraw[0] + 1 )
			rename w, condtemp
			wave condtemp
				
		elseif( i == scanbinsraw[0] + 2 )
			rename w, sampsdev
			wave sampsdev
			
		elseif( i == scanbinsraw[0] + 3 )
			rename w, shsdev
			wave shsdev
					
		elseif( i == scanbinsraw[0] + 4 )	
			rename w, cpcflowflag
			wave cpcflowflag
				
		elseif( i == scanbinsraw[0] + 5 )
			rename w, cpcfillflag
			wave cpcfillflag
			
		elseif( i == scanbinsraw[0] + 6 )
			rename w, shflowflag
			wave shflowflag
				
		elseif( i == scanbinsraw[0] + 7 )
			rename w, numflagraw
			wave numflagraw
														
		endif
	endfor
	
	for( i = 0; i < 8; i += 1 )		//rename bin count autonamed waves
//	for( i = 0; i < scanbinsraw[0]; i += 1 )		//rename bin count autonamed waves

		wavenum = 2*scanbinsraw[0] + 13 + i
//		wavestr = "wave" + num2str(wavenum)
//		wave w0 = $wavestr
		wavestr1 = "wave" + num2str(wavenum)// + scanbinsraw[0])
		wave w = $wavestr1
		
		if( i == 0 )
			rename w, sattemp
			wave sattemp
			
		elseif( i == 1 )
			rename w, condtemp
			wave condtemp
				
		elseif( i == 2 )
			rename w, sampsdev
			wave sampsdev
			
		elseif( i == 3 )
			rename w, shsdev
			wave shsdev
					
		elseif( i == 4 )	
			rename w, cpcflowflag
			wave cpcflowflag
				
		elseif( i == 5 )
			rename w, cpcfillflag
			wave cpcfillflag
			
		elseif( i == 6 )
			rename w, shflowflag
			wave shflowflag
				
		elseif( i == 7 )
			rename w, numflagraw
			wave numflagraw
														
		endif
	endfor
	
	for( i = 0; i < numscans ; i += 1 )				//reorg bincnt[scan] waves into scan[bin] waves
		scanstr = "nm_scan" + num2str( i + 1)
		scanstr1 = "cnt_scan" + num2str( i + 1)
		make /o/n =(scanbinsraw[i]) $scanstr
		make /o/n =(scanbinsraw[i]) $scanstr1
		wave scandia = $scanstr
		wave scancnt = $scanstr1
				
		for( j = 0; j < scanbinsraw[i]; j += 1 )
		
			binstr1 = "bindia" + num2str( j + 1 )
			binstr2 = "bincnt" + num2str( j + 1 )
			wave bindia = $binstr1
			wave bincnt = $binstr2
			scandia[j] = bindia[i]
			scancnt[j] = bincnt[i]
			
		endfor
	endfor
	
End


//**************************************************************************************************
// assumes particular timing
// Ambient = starts on zeros
// TD = starts on 5's
Function SplitSEMS_byTDstate()
	
//	setdatafolder root:smps:Merged
	string dt_smps_str="dt_smps"
	string dndlogdp_str="dNdlogDp_smps"
	variable tref_ambient=0
	variable tref_denuded=2
	variable tref_treatspecial1=25
	variable tref_treatspecial2=55
	prompt dt_smps_str, "Select time wave", popup WaveList("*",";","")
	prompt dndlogdp_str, "Select number distribution wave",popup WaveList("*",";","")
	prompt tref_ambient, "Starting minute of ambient data"
	prompt tref_denuded, "Starting minute of TD data"
	prompt tref_treatspecial1, "Minute where only 1/2 a scan was completed"
	prompt tref_treatspecial2, "Minute where only 1/2 a scan was completed"
	DoPrompt "Select waves", dt_smps_str, dndlogdp_str,tref_ambient,tref_denuded//, tref_treatspecial1, tref_treatspecial2
	
	wave dt_SMPS = $dt_smps_str
	wave dndlogdp_SMPS = $dndlogdp_str
	
	if(numpnts(dt_smps)!=dimsize(dndlogdp_smps,0))
		Abort "Waves must have same length. Aborting."
	endif
	
	variable V_time = numpnts(dt_SMPS)
	make/o/d/n=(V_time) maskTD = nan, maskSpecial=1
	variable i
	string dt_string, dt, dtt
	
	for(i=0;i<V_time;i+=1)
		dt_string = secs2time(dt_smps[i],2)
		dt = dt_string[4]//stringfromlist(1,dt_string,":")
		dtt = dt_string[3,4]
		if(str2num(dt)==tref_ambient)
			maskTD[i] = 0
		elseif(str2num(dt)==tref_denuded)
			maskTD[i] = 1
		endif
//		if(str2num(dtt)==tref_treatspecial1 || str2num(dtt)==tref_treatspecial2)
//			maskSpecial[i] = 2
//		else
//			maskSpecial[i] = 1
//		endif
	endfor
	
	string TD = dndlogdp_str + "_TD"
	string RT = dndlogdp_str + "_RT"
	duplicate/o dndlogdp_SMPS $TD, $RT
	wave dndlogdp_TD = $TD
	wave dndlogdp_RT = $RT
	dNdlogDp_TD[][] = maskTD[p]==1 ? dndlogdp_SMPS[p][q] : nan
	dNdlogDp_RT[][] = maskTD[p]==0 ? dndlogdp_SMPS[p][q] : nan
//	dNdlogDp_TD[][] *= maskSpecial[p]
	
//	killwaves/z maskTD, maskSpecial
end
//////////////////////////////////////////// SEMS ABOVE ///////////////////////////////
