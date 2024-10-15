/* dem2tga.c (c) 2004 Michael Hogsett
 *  derived from :
 * 
 * dem2tga.c (c) 1997 Jon Larimer
 *
 *  This is a program to convert USGS Digital Elevation Model (DEM) data
 *  into a .TGA file, viewable by many graphics applications.  This program
 *  is public domain, do whatever you want with it, except claim it as your
 *  own or sell it. If you find a way to make money from my work, I want a
 *  piece of the action.
 *
 *  This should compile with any ANSI standard compiler, including Unix
 *  and DOS versions of gcc and many other C compilers.  If you have any
 *  questions or comments email jonl@alltel.net.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/errno.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#define kTYPE_A_SIZE 1024
#define kTYPE_B_SIZE 8192

#define kINT_LENGTH 6
#define kFLOAT_LENGTH 12
#define kDOUBLE_LENGTH 24

void usage() { 
  fprintf(stderr,"usage: dem2tga [-v] [-n] -m min_elev -s scale_factor dem_file tga_file\n");
  fprintf(stderr,"       dem2tga [-v] [-n] -e dem_file1 [dem_file2 ...]\n");
  fprintf(stderr,"                -v : verbose\n");
  fprintf(stderr,"                -n : print DEM name(s) to stdout\n");
  fprintf(stderr,"                -e : extract global elevation scale\n");
  fprintf(stderr,"                -m n -s n : force min_elev to n and scale to n\n");
  exit(1);
}

int writetgaheader(FILE *fptr, int r, int c) {
  int i;
  fputc(0, fptr); fputc(1, fptr); fputc(1, fptr);
  fputc(0, fptr); fputc(0, fptr); fputc(0, fptr);
  fputc(1, fptr); fputc(24, fptr); fputc(0, fptr);
  fputc(0, fptr); fputc(0, fptr); fputc(0, fptr);
  fputc((char)(c & 0x00ff), fptr); 
  fputc((char)((c & 0xff00) >> 8), fptr);
  fputc((char)(r & 0x00ff), fptr);
  fputc((char)((r & 0xff00) >> 8), fptr);
  fputc(8, fptr); fputc(0, fptr);

  for(i=0; i<=255; i++) { 
	fputc(i, fptr); 
	fputc(i, fptr);
	fputc(i, fptr); 
  }
  fflush(fptr);
  return 0;
}

int getnextint(char *s) { 
  char buf[kINT_LENGTH+1];
  int i;
  for(i=0;i<kINT_LENGTH;i++) { buf[i] = s[i]; }
  buf[kINT_LENGTH] = '\0';
  return(atoi(buf));
}

double getnextdouble(char *s) { 
  char buf[kDOUBLE_LENGTH+1];
  int i;
  for(i=0;i<kDOUBLE_LENGTH;i++) { 
	buf[i] = s[i]; 
	if( buf[i] == 'D' ) { 
	  buf[i] = 'E';
	}
  }
  buf[kDOUBLE_LENGTH] = '\0';
  return(atof(buf));
}

float getnextfloat(char *s) { 
  char buf[kFLOAT_LENGTH+1];
  int i;
  for(i=0;i<kFLOAT_LENGTH;i++) { 
	buf[i] = s[i]; 
	if( buf[i] == 'D' ) { 
	  buf[i] = 'E';
	}
  }
  buf[kFLOAT_LENGTH+1] = '\0';
  return((float)atof(buf));
}

int main(int argc, char **argv) {
  FILE *demfile;
  FILE *tgafile;
  char name[145], type_a_record[1025], type_b_record[8193];
  int dem_level_code, pattern_code, plan_ref_sys_code, zone_code, accuracy_code;  
  int ground_units_code, elev_units_code,poly_sides, profile_dim, profile_num;
  int profile_elevs, tga_dim_x, tga_dim_y, current_profile, this_profile_dim, elev;
  int this_profile_id, this_profile_elevs, this_profile_columns, i, offset;
  float x_res, y_res, z_res;
  double map_proj_param[15], poly_verts[8], width, height, this_profile_long, this_profile_lat;
  double min_elev, max_elev, elev_range, scaling_factor, angle_from_axis, this_profile_local_elev;
  double this_profile_min_elev, this_profile_max_elev;

  int verbose, ch, dump_header, scale_provided, min_elev_provided, output_location;
  double scale, elev_extract, provided_elev;
  double global_min_elev, global_max_elev;
  int first_min_elev = 0;
  global_min_elev = 0.0;
  global_max_elev = 0.0;
  output_location = 0;

  verbose = 0;
  dump_header = 0;
  scale_provided = 0;
  min_elev_provided = 0;
  elev_extract = 0;

  while ((ch = getopt(argc, argv, "delm:ns:v")) != -1)
	switch(ch) { 
	case 'e':
	  elev_extract=1;
	  if ( verbose == 1 ) { fprintf(stderr,"Extracting Elevations from given files\n"); }
	  break;
	case 'l':
	  output_location = 1;
	  if ( verbose == 1 ) { fprintf(stderr,"Writing DEM Location to stdout\n"); }
	  break;	  
	case 'm':
	  if ( min_elev_provided != 0 ) { 
		fprintf(stderr,"Error: min elev already provided\n");
		exit(1);
	  }
	  provided_elev = atof(optarg);
	  if ( verbose == 1 ) { fprintf(stderr,"Using given min elev %.8f\n",provided_elev); }
	  min_elev_provided = 1;
	  break;
	case 'n':
	case 'd':
	  dump_header=1;
	  if ( verbose == 1 ) { fprintf(stderr,"Writing DEM name to stdout\n"); }
	  break;
	case 's':
	  if ( scale_provided != 0 ) { 
		fprintf(stderr,"Error: scale already provided\n");
		exit(1);
	  }
	  if (( scale = atof(optarg)) <= 0.0 ) {
		fprintf(stderr,"Error : scale less than or equal to zero. \"%s\"\n",optarg);
		exit(1);
	  } else { 
		if ( verbose == 1 ) { fprintf(stderr,"Using given scaling factor %.8f\n",scale); }
		scale_provided = 1;
	  }
	  break;
	case 'v':
	  verbose=1;
	  if ( verbose == 1 ) { fprintf(stderr,"Verbose set\n"); }
	  break;
	case '?':
	default:
	  usage();
	}
  argc -= optind;
  argv += optind;

  // expect an input and an output file.
  if ( elev_extract == 0 && dump_header != 1 ) {
	if(argc < 2) {
	  usage();
	}
  } else { 
	if (argc < 1) { 
	  usage();
	}
  }

  if ( ((scale_provided == 1) ^ (min_elev_provided == 1)) ) { 
	fprintf(stderr,"Error:  need both scale and min_elev.  Exiting.\n");
	exit(1);
  }

  if ( elev_extract == 1) { 
	/* for each argv, extract elevation extremes, update global
	   elevations, then calculate scaling factor based on elevation
	   extremes and output scaling factor
	*/
	i=0;
	while ( i < argc ) { 
	  if ( verbose == 1 ) { fprintf(stderr,"DEM File %s: ",argv[i]); }
	  if((demfile = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "Error : %s.  Exiting.\n",strerror(errno));
		exit(1);
	  }
	  fread(type_a_record,sizeof(char),kTYPE_A_SIZE,demfile);
	  fclose(demfile);
	  type_a_record[kTYPE_A_SIZE] = '\0';
	  min_elev = getnextdouble(&type_a_record[738]);
	  max_elev = getnextdouble(&type_a_record[762]);
	  elev_range = max_elev - min_elev;
	  if ( verbose == 1 ) { 
		fprintf(stderr, " min,max elev: %.2f to %.2f, range %.2f\n", min_elev,max_elev,elev_range);
	  }
	  if ( elev_range < 0.0 ) { 
		fprintf(stderr,"Negative elevation range.  Exiting.\n");
		exit(1);
	  }
	  if ( first_min_elev == 0 ) { 
		global_min_elev = min_elev;
		first_min_elev = 1;
	  } else {
		if ( min_elev < global_min_elev ) { 
		  global_min_elev = min_elev;
		}
	  }
	  if ( max_elev > global_max_elev ) { 
		global_max_elev = max_elev;
	  }
	  i++;
	}
	elev_range = global_max_elev - global_min_elev;
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Global min,max elev: %.2f to %.2f, range %.2f\n", global_min_elev,global_max_elev,elev_range);
	}
	if ( elev_range == 0.0 ) { 
	  scaling_factor = 0.0;
	}
	else { 
	  scaling_factor = 255.0 / elev_range;
	}
	fprintf(stdout,"%g %g\n",global_min_elev,scaling_factor);
	exit(0);
  } else if ( dump_header == 1 ) { 
	// dump name
	// dump last lat,long pair
	i=0;
	while ( i < argc ) { 
	  int j=0;
	  if((demfile = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "Error : %s.  Exiting.\n",strerror(errno));
		exit(1);
	  }
	  fread(type_a_record,sizeof(char),kTYPE_A_SIZE,demfile);
	  fclose(demfile);
	  type_a_record[kTYPE_A_SIZE] = '\0';
	  for(j=0; j <= 39; j++) { name[j] = type_a_record[j]; }
	  name[40] = '\0';
	  for(j=39; j > 0; j--) {
		if(!isspace(name[j])) { j=0; } else { name[j] = '\0'; }
	  }
	  fprintf(stdout,"NAME=\"%s\";\n",name);
	  for(i=0; i < 8; i++) {
		poly_verts[i] = getnextdouble(&type_a_record[546 + i * kDOUBLE_LENGTH]);
	  }
	  if ( poly_verts[7] < 0.0 ) { 
		fprintf(stdout,"LOC=\"%.2fS,",fabs(poly_verts[7]/3600));
	  } else { 
		fprintf(stdout,"LOC=\"%.2fN,",fabs(poly_verts[7]/3600));
	  }
	  if ( poly_verts[6] < 0.0 ) { 
		fprintf(stdout,"%.2fW\";\n",fabs(poly_verts[6]/3600));
	  } else { 
		fprintf(stdout,"%.2fE\";\n",fabs(poly_verts[6]/3600));		
	  }
	  i++;
	}
	exit(0);
  } else  { 

	// fprintf(stderr,"%s %s\n",argv[0],argv[1]);
	if((demfile = fopen(argv[0], "r")) == NULL) {
	  fprintf(stderr, "Error : %s.  Exiting.\n",strerror(errno));
	  exit(1);
	}
	
	/* DEM Type A Records */
	
	// slurp in the Type A Record header.  It consists of the first 1024
	// bytes of the DEM file.
	
	fread(type_a_record,sizeof(char),kTYPE_A_SIZE,demfile);
	type_a_record[kTYPE_A_SIZE] = '\0';
	
	//  fprintf(stderr,"### RAW HEADER ####\n");
	//  fprintf(stderr,"%s\n",type_a_record);
	//  fprintf(stderr,"### END HEADER ####\n");
	
	/*
	  All fields are represented in ASCII
	  chars are the literal ascii string
	  flags are 2 bytes 
	  shorts are 4 bytes 
	  ints are 6 bytes 
	  floats are 12 bytes
	  doubles are 24 bytes
	*/
	
	/* Field 1 - char string DEM name field.  bytes 0 to 143 
	   used only if outputting name
	*/
	if ( dump_header == 1 ) { 
	  char cur_c, prev_c;
	  for(i=0; i <= 39; i++) { name[i] = type_a_record[i]; }
	  name[40] = '\0';
	  for(i=39; i > 0; i--) {
		if(!isspace(name[i])) { i=0; } else { name[i] = '\0'; }
	  }
	  
	  // REMOVE all ocurrences of 2 spaces in a row
	  cur_c = '.'; 
	  prev_c = '.';
	  i = 0;
	  while ( cur_c != '\0' && i <= 39 ) { 
		cur_c = name[i];
		if ( cur_c == ' ' && prev_c == ' ' ) { 
		  bcopy(&name[i],&name[i-1],strlen(&name[i]));
		  name[strlen(name)-1] = '\0';
		  i=0;
		  prev_c = '.';
		} else { 
		  prev_c = cur_c;
		  i++;
		}
	  }
	  // REPLACE ' - ' with '-'
	  i = 0;
	  while ( name[i] != '\0' && i < (strlen(name) - 2) ) { 
		if ( name[i] == ' ' && name[i+1] == '-' && name[i+2] == ' ' ) {
		  // replace name[i] with '-'
		  name[i] = '-';
		  // move name[i+3] to name[i+1]
		  // add nulls to end
		  bcopy(&name[i+3],&name[i+1],strlen(&name[i+3]));
		  name[strlen(name)-2] = '\0';
		  //reset i
		  i = 0;
		} else { 
		  i++;
		}
	  }
	  // REPLAE ' ' with '_'
	  i = 0;
	  while ( name[i] != '\0' && i < (strlen(name)) ) { 
		if ( name[i] == ' ' ) { 
		  name[i] = '_';
		}
		i++;
	  }
	  fprintf(stdout, "DEM Name:\"%s\"\n", name);
	}
	
	/* Field 2 - int DEM level code bytes 144 to 149 
	   unused
	*/
	if ( verbose == 1 ) { 
	  dem_level_code = getnextint(&type_a_record[144]);
	  fprintf(stderr, "DEM Level Code %d ",dem_level_code);
	  if ( dem_level_code == 3 ) { 
		fprintf(stderr,"(processed by DMA)\n");
	  } else { 
		fprintf(stderr,"(unknown)\n");
	  }
	}

	/* Field 3 - int pattern code bytes 150 to 155
	   unused
	*/
	if ( verbose == 1 ) { 
	  pattern_code = getnextint(&type_a_record[150]);
	  fprintf(stderr, "Pattern Code %d ",pattern_code);
	  if ( pattern_code == 1 ) { 
		fprintf(stderr,"(regular elevation pattern)\n");
	  } else { 
		fprintf(stderr,"(unknown)\n");
	  }
	}

	/* Field 4 - int planimetric reference system code bytes 156 to 161
	   unused
	*/
	if ( verbose == 1 ) {   
	  plan_ref_sys_code = getnextint(&type_a_record[156]);
	  fprintf(stderr, "Planimetric Reference System Code %d ",plan_ref_sys_code);
	  if ( plan_ref_sys_code == 0 ) { 
		fprintf(stderr,"(geographic coordinate system)\n");
	  } else { 
		fprintf(stderr,"(unknown)\n");
	  }
	}

	/* Field 5 - int zone code bytes 162 to 167
	   unused
	*/
	if ( verbose == 1 ) { 
	  zone_code = getnextint(&type_a_record[162]);
	  fprintf(stderr, "Zone Code %d ",zone_code);
	  if ( zone_code != 0 ) { 
		fprintf(stderr,"- Unexpected Value");
	  }
	  fprintf(stderr,"\n");
	}

	/* Field 6 - double x 15 map projection parameters bytes 168 to 527
	   unused
	*/
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Map projection parameters ");
	  for(i=0; i < 15; i++) { 
		map_proj_param[i] = getnextdouble(&type_a_record[168+i*kDOUBLE_LENGTH]);
		fprintf(stderr,".");
	  }
	  fprintf(stderr,"\n");
	}
	
	/* Field 7 - int ground units code bytes 528 to 533
	 */
	ground_units_code = getnextint(&type_a_record[528]);
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Ground Units Code %d ",ground_units_code);
	  if ( ground_units_code == 3 ) { 
		fprintf(stderr, "(arc-seconds)\n");
	  } else { 
		fprintf(stderr, "(unknown)\n");
	  }
	}

	/* Field 8 - int elevation units code bytes 534 to 539
	 */
	elev_units_code = getnextint(&type_a_record[534]);
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Elevation Units %i ", elev_units_code);
	  if ( elev_units_code == 2 ) {
		fprintf(stderr, "(meters)\n");
	  } else { 
		fprintf(stderr, "(unknown)\n");
	  }
	}

	/* Field 9 - int dem polygon sides bytes 540 545
	 */
	poly_sides = getnextint(&type_a_record[540]);
	if ( verbose == 1 ) { 
	  fprintf(stderr,"Number of sides of the DEM polygon: %d\n",poly_sides);
	}
	if ( poly_sides != 4 ) {
	  fprintf(stderr,"I don't know how to deal with non-rectangluar DEMs.  Exiting.\n");
	  exit(1);
	} 
  
	/* Field 10 - (double,double) x 4 polygon vertex coords bytes 546 to 737
	 */
	if ( verbose == 1 ) {   fprintf(stderr,"Ground coordinates of 4 corners of DEM: "); }
	for(i=0; i < 8; i++) {
	  poly_verts[i] = getnextdouble(&type_a_record[546 + i * kDOUBLE_LENGTH]);
	  if ( verbose == 1 ) { 
		if ( (i%2) == 0 ) {
		  fprintf(stderr,"(%.0f,",poly_verts[i]/3600.0);
		} else { 
		  fprintf(stderr,"%.0f) ",poly_verts[i]/3600.0);
		}
	  }
	}
	width = fabs(poly_verts[0] - poly_verts[4]);
	height = fabs(poly_verts[1] - poly_verts[3]);

	if ( verbose == 1 ) { 
	  fprintf(stderr,"\nDEM ground coordinate width,height : %.2f,%.2f\n",width,height);
	}
	
	if ( output_location == 1 ) { 
	  fprintf(stdout,"DEM Location:(%.4f,%.4f)-",poly_verts[0]/3600,poly_verts[1]/3600);
	  fprintf(stdout,"(%.4f,%.4f)\n",poly_verts[4]/3600,poly_verts[5]/3600);
	}

	/* Field 11 double,double min and max elevations in DEM bytes 738 to
	   761 and 762 to 785
	*/
	min_elev = getnextdouble(&type_a_record[738]);
	max_elev = getnextdouble(&type_a_record[762]);
	elev_range = max_elev - min_elev;
	if ( verbose == 1 ) { 
	  fprintf(stderr, "DEM min,max elevation: %.2f to %.2f, range %.2f\n", min_elev,max_elev,elev_range);
	}
	if ( elev_range < 0.0 ) { 
	  fprintf(stderr,"Negative elevation range.  Exiting.\n");
	  exit(1);
	}
	if ( min_elev_provided == 1 ) {
	  min_elev = provided_elev;
	  if ( verbose == 1 ) { 
		fprintf(stderr, "Using Provided min elev : %.8f\n", min_elev);
	  }
	}
  
	/* this is the value to multiply elevations by to get a number
	   between 0 and 255 to determine the color of the pixel */
	if ( elev_range == 0.0 ) { 
	  scaling_factor = 0.0;
	} else { 
	  scaling_factor = 255.0 / elev_range;
	}
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Calculated Scaling factor : %.5f\n", scaling_factor);
	}
	if ( scale_provided == 1 ) { 
	  scaling_factor = scale; 
	  if ( verbose == 1 ) { 
		fprintf(stderr, "Using Provided Scaling factor : %.5f\n", scaling_factor);
	  }
	} 
	
	/* Field 12 double ccw angle from primary axis bytes 786 to 809
	   unused
	*/
	angle_from_axis = getnextdouble(&type_a_record[786]);
	if ( verbose == 1 ) { 
	  fprintf(stderr, "CCW Angle from primary axis: %.2f\n",angle_from_axis);
	}

	/* Field 13 int accuracy code bytes 810 to 815
	   unused
	*/
	accuracy_code = getnextint(&type_a_record[810]);
	if ( verbose == 1 ) { 
	  fprintf(stderr, "Accuracy code %d, ",accuracy_code);
	  if ( accuracy_code == 0 ) { 
		fprintf(stderr,"no class C records follow.\n");
	  } else {
		fprintf(stderr,"class C records follow.\n");
	  }
	}

	/* Field 14 float x 3 DEM spatial resolution bytes 816 to 851
	 */
	x_res = getnextfloat(&type_a_record[816]);
	y_res = getnextfloat(&type_a_record[828]);
	z_res = getnextfloat(&type_a_record[840]);
	if ( verbose == 1 ) { 
	  fprintf(stderr,"X, Y, Z Spatial Resolution : %.2f, %.2f, %.2f ",x_res,y_res,z_res);
	  if ( ground_units_code == 3 ) { 
		fprintf(stderr,"arcseconds");
	  } 
	  fprintf(stderr,"\nExpecting %d profiles ",(int)(height/x_res+1));
	  fprintf(stderr,"of %d elevations each.\n",(int)(width/y_res+1));
	}
	
	/* Field 15 int x 2 profile array rows and columns bytes 852 to 857
	   and 858 to 863
	*/
	profile_dim = getnextint(&type_a_record[852]);
	profile_num = getnextint(&type_a_record[858]);
	if ( verbose == 1 ) { 
	  fprintf(stderr,"There are %d %d dimensional profiles in this DEM.\n",profile_num,profile_dim);
	}
	if ( profile_dim != 1 ) { 
	  fprintf(stderr, "Can't handle %d dimensional DEM files.  Exiting\n",profile_dim);
	  exit(1);
	}
	if ( profile_num != (int)(height/x_res+1) ) { 
	  fprintf(stderr, "Unexpected number of profiles.  Exiting.\n");
	  exit(1);
	}
	
	if ( verbose == 1 ) { 
	  fprintf(stderr,"Discarding remaining Type A Record fields.\n");
	}
	
	// dive into first Type B Record to retrieve
	// the rows(elevations) per profile
	if ( verbose == 1 ) { 
	  fprintf(stderr,"Entering First Type B Record to get elevations per profile: ");
	}
	fseek(demfile,1036,SEEK_SET);
	fread(type_a_record,sizeof(char),6,demfile);
	// RE-USE type_a_record buffer space
	type_a_record[6] = '\0';
	profile_elevs = getnextint(&type_a_record[0]);
	if ( verbose == 1 ) { fprintf(stderr," %d\n",profile_elevs); }
	
	tga_dim_x = profile_elevs;
	tga_dim_y = profile_num;
	
	if ( verbose == 1 ) { 
	  fprintf(stderr,"TGA Image is %d x %d \n",tga_dim_x,tga_dim_y);
	  fprintf(stderr,"Writing TGA header\n");
	}
	if((tgafile = fopen(argv[1], "wb+")) == NULL ) { 
	  fprintf(stderr, "%s: fopen: %s", argv[0], strerror(errno));
	  exit(1);
	}
	writetgaheader(tgafile, tga_dim_y, tga_dim_x);
	
	if ( verbose == 1 ) { 
	  fprintf(stderr,"Setting File Offset to first Type B Record\n");
	}
	fseek(demfile,kTYPE_A_SIZE,SEEK_SET);
	fflush(stderr);

	/***************************************************************************** 
	 * DEM Type B Records  
	 *****************************************************************************/
	if ( verbose == 1 ) { 
	  fprintf(stderr,"Parsing DEM Type B Records (dot every 100 profiles): ");
	}
	current_profile = 1;
	while( current_profile <= profile_num ) {
	  i = current_profile -1;
	  if ( verbose == 1 ) { 
		if ( (i > 0) && ((i%100) == 0) ) { 
		  fprintf(stderr,".");
		}
	  }

	  // read type b record into 8k buffer
	  fread(type_b_record,sizeof(char),kTYPE_B_SIZE,demfile);
	  type_b_record[kTYPE_B_SIZE] = '\0';
	  
	  /* Field 1 int x 2
		 profile row and column id
		 column appears to be profile id
		 row appears to be this profiles dimension number
		 bytes 0 to 5 and 6 to 11
	  */
	  this_profile_dim =  getnextint(&type_b_record[0]);
	  this_profile_id = getnextint(&type_b_record[6]);
	  if ( this_profile_id != current_profile ) { 
		fprintf(stderr,"Expecting profile id %d, got %d.  Exiting.\n",current_profile,this_profile_id);
		exit(1);
	  }
	  if ( this_profile_dim != profile_dim ) { 
		fprintf(stderr,"Expecting profile dimension %d, got %d.  Exiting.\n",profile_dim,this_profile_dim);
		exit(1);
	  }
	  /* Field 2 int x 2 
		 profile rows and columns
		 rows is elevations (samples) in this profile
		 columns  matches dimensions?	   
		 bytes 12 to 17 and 18 to 23
	  */
	  this_profile_elevs = getnextint(&type_b_record[12]);
	  this_profile_columns = getnextint(&type_b_record[18]);
	  if ( this_profile_elevs != profile_elevs ) { 
		fprintf(stderr,"Expecting %d elevations, got %d.  Exiting.\n",this_profile_elevs,profile_elevs);
		exit(1);
	  }
	  if ( this_profile_columns != 1 ) { 
		fprintf(stderr,"Expecting 1 profile columns, got %d.  Exiting.\n",this_profile_columns);
		exit(1);
	  }
	  
	  /* Field 3 double x 2
		 ground coords of first elevation in profile
		 bytes 24 to 47 and 48 to 71 
	  */
	  this_profile_long = getnextdouble(&type_b_record[24]);
	  this_profile_lat = getnextdouble(&type_b_record[48]);
	  /* !!! check to make sure these are not out of bounds of the information in the DEM header */
	  
	  /* Field 4 double 
		 profile local datum elevation
		 always 0.0 (sealevel) for 1 degree DEM
		 bytes 72 to 95
	  */
	  this_profile_local_elev = getnextdouble(&type_b_record[72]);
	  if ( this_profile_local_elev != 0.0 ) { 
		fprintf(stderr,"Expecting 0.0, got %.1f.  Exiting.\n",this_profile_local_elev);
		exit(1);
	  } 
	  
	  /* Field 5 double x 2 
		 min and max elevations for this profile 
		 bytes 96 to 119 and 120 to 143
	  */
	  this_profile_min_elev = getnextdouble(&type_b_record[96]);
	  this_profile_max_elev = getnextdouble(&type_b_record[120]);
	  /* !!! check to make sure these are not out of bounds of the information in the DEM header */
	  
	  /* Field 6 int x profile_elevs
		 elevation samples
		 byte 144 to the end
	  */
	  offset = 144;
	  for(i=0; i < profile_elevs; i++) { 
		// correct at 1k boundaries
		if ( i>0 && ((i-146) % 170 == 0) ) {
		  offset += 4;
		}
		elev = getnextint(&type_b_record[offset+(i*kINT_LENGTH)]);
		fputc((int)( (elev - min_elev) * scaling_factor),tgafile);
	  }
	  current_profile++;
	}
	if ( verbose == 1 ) { fprintf(stderr, " done.\n"); }
	fclose(tgafile);
	fclose(demfile);
	return(0);
  }
  return(0);
}
