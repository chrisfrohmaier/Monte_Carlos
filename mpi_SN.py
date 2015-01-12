#Test
from mpi4py import MPI
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from matplotlib.patches import Circle
from shapely.geometry import Polygon, Point
from descartes.patch import PolygonPatch
from astropy.time import Time
import subprocess, brewer2mpl, math, sncosmo, psycopg2
colors = brewer2mpl.get_map('Set3', 'Qualitative', 12).mpl_colors
import time as TI

##Connecting to the Database
conn = psycopg2.connect(host='srv01050.soton.ac.uk', user='frohmaier', password='rates', database='frohmaier')
cur = conn.cursor()
conn2 = psycopg2.connect(host='srv01050.soton.ac.uk', user='frohmaier', password='rates', database='frohmaier')
cur2=conn2.cursor()
conn3 = psycopg2.connect(host='srv01050.soton.ac.uk', user='frohmaier', password='rates', database='frohmaier')
cur3 = conn3.cursor()
cur.execute("SELECT DISTINCT ON (ujd) ujd from obs ORDER BY ujd;")
d=cur.fetchall()
dates=[]
for i in range(len(d)):
	dates.append([int(i), float(d[i][0])])
all_data=np.array(dates)
#print all_data	
cur.close()
def Create_Date_Array(peak_date,date_array):
	date_a, index=Find_Nearest_Date(peak_date,date_array)
	
	Sn_Dates=[]
	min_array, min_date_index=Find_Nearest_Date(peak_date-19,date_array)
	max_array, max_date_index=Find_Nearest_Date(peak_date+51, date_array)
	for i in range(min_date_index,max_date_index):
		#print i
		Sn_Dates.append(date_array[i][1])
	return Sn_Dates
def Find_Nearest_Date(date_ujd,nupy_array):
	idx=(np.abs(nupy_array[:,1]-date_ujd)).argmin()
	#print 'Index: ', idx
	return nupy_array[idx], idx
def Random_Gen_RA_DEC_MJD(date):
	Good=False
	while Good==False:
		RA=np.random.uniform(-0.52,-0.452)
		DEC=np.random.uniform(-0.350,-0.278)
		cord=Point(RA,DEC)
		area_p=Polygon([(-0.482, -0.350),(-0.453, -0.331),(-0.48722712, -0.2789),(-0.516, -0.2972)])
		if cord.within(area_p)==True:
			Peak_Date=np.random.uniform(min(date[:,1])-50,max(date[:,1])+20)
			Good=True
			return RA, DEC, Peak_Date
def pol_2_cart(ra,dec):
	ra_r=np.radians(ra)
	dec_r=np.radians(dec)
	x=np.cos(dec_r)*np.cos(ra_r)#*math.sin(math.pi/2.0-dec_r)
	#y=math.sin(ra_r)*math.sin(math.pi/2.0-dec_r)
	y=np.cos(dec_r)*np.sin(ra_r)
	#print x,y
	return x,y
def cart_2_pol(x,y):
	ra_r=np.arctan(y/x)
	dec_r=np.arccos(np.sqrt(y**2. + x**2.))

	ra=np.degrees(ra_r)
	dec=np.degrees(dec_r)

	return ra+180., dec
def which_ccd(ujd, ra, dec):
	cur2.execute("SELECT  * from obs where ujd=%s;",(float(ujd),))
	m=cur2.fetchall()
	for ln in m:
		ujd=float(ln[0])
		seeing_new=float(ln[1])
		ub1_zp_new=float(ln[2])
		lmt_mag_new=float(ln[3])

		ccdid=int(ln[4])
		goodpixarea=float(ln[5])
		ra_ul=float(ln[6])
		dec_ul=float(ln[7])
		x_ul,y_ul=pol_2_cart(ra_ul,dec_ul)
		ra_ur=float(ln[8])
		dec_ur=float(ln[9])
		x_ur,y_ur=pol_2_cart(ra_ur,dec_ur)
		ra_lr=float(ln[10])
		dec_lr=float(ln[11])
		x_lr,y_lr=pol_2_cart(ra_lr,dec_lr)
		ra_ll=float(ln[12])
		dec_ll=float(ln[13])
		x_ll,y_ll=pol_2_cart(ra_ll,dec_ll)
		
		ccd_polygon=Polygon([(x_ul,y_ul),(x_ur,y_ur),(x_lr,y_lr),(x_ll,y_ll)])

		sn_object= Point((ra, dec))

		if ccd_polygon.contains(sn_object)==True:
			return ccdid, 1-(goodpixarea/0.6603)
	return	99.9, 1.

bpass=np.loadtxt('PTF48R.dat')
wavelength=bpass[:,0]
transmission=bpass[:,1]
band=sncosmo.Bandpass(wavelength,transmission, name='ptf48r')
sncosmo.registry.register(band, force=True)

def Gen_SN(peak_time, redshift, colour,x_1, date):
	source=sncosmo.get_source('salt2',version='2.4') #Importing SALT2 Model
	model=sncosmo.Model(source=source) 

	alpha=0.141
	beta=3.101
	int_dis=np.random.normal(0.,0.15)

	mabs= -19.05 - alpha*x_1 + beta*colour + int_dis
	model.set(z=redshift,t0=peak_time,x1=x_1, c=colour) #Setting redshift
	
	model.set_source_peakabsmag(mabs,'bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
	#model.set(x1=x_1, c=colour)

	band=sncosmo.get_bandpass('ptf48r') #Retrieving the ptf48r bandpass 
    
	time=Create_Date_Array(peak_time,date) #setting an arbitrary time span to cover the model
	
	maglc=model.bandmag('ptf48r','ab',time) #Creating a magnitude array of the lightcurve  
	fluxlc=model.bandflux('ptf48r',time) #Creating a flux array of the lightcurve
	absmagb=model.source_peakabsmag('bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
	absmag_r=model.source_peakabsmag('ptf48r','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
	return peak_time, time, maglc, fluxlc, absmagb, absmag_r, redshift, colour, x_1, int_dis #returns

nproc = MPI.COMM_WORLD.Get_size()   
# number of processes
my_rank = MPI.COMM_WORLD.Get_rank()   
# The number/rank of this process
my_node = MPI.Get_processor_name()    
# Node where this MPI process runs

N_MODELS_TOTAL = 1000  # total number of models to run
n_models = N_MODELS_TOTAL / nproc
# number of models for each thread

my_nmin = (my_rank * n_models) 
my_nmax = my_nmin + n_models

time_init = TI.time()
for i in range( my_nmin, my_nmax): 
	ra,dec,peak_date=Random_Gen_RA_DEC_MJD(all_data)
	xone=np.random.uniform(-3.0,3.0)
	color=np.random.uniform(-0.3,0.3)
	zedshift=np.random.uniform(0.0,0.1)

	ra1,dec1=cart_2_pol(ra,dec)
	#print zedshift

	peak_time, time, maglc, fluxlc, absmagb, absmag_r, redshift, colour, x_1, int_dis=Gen_SN(peak_date, zedshift, color,xone, all_data)

##Inserting the Values into the sne table 

	#print 'Doing CCD Shit'
	probs=[]
	t_ccd=[]
	for j in range(0,len(time)):
		#print time[i]
		#this_ccd, not_detect_prob=which_ccd(time[j], ra, dec)
		probs.append(0)
		t_ccd.append(this_ccd)
		#print this_ccd
		#print probs	
	total_prob_not=np.prod(probs)
	cur = conn.cursor()

	#print maglc, time
	cur.execute("INSERT INTO sne (peak_date, ra, dec, absolute_mag, redshift, x1, color, int_dispersion, prob_not) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s) returning sn_id",(peak_time,ra1,dec1,absmag_r,redshift,x_1, colour, int_dis,total_prob_not))
	sn_id=cur.fetchone()
	print 'This is the SN ID: ', sn_id[0]
	for k in range(0, len(time)):
		cur3.execute("INSERT INTO sim_epoch (sn_id, ujd, ra, dec, magnitude, ccd, prob_not_detected) VALUES (%s,%s,%s,%s,%s,%s,%s)", (sn_id[0], float(time[k]), ra1, dec1, maglc[k], t_ccd[k], probs[k]))
		conn3.commit()
	
	#print len(time), len(probs)
	#print 'Total Not Prob: ', total_prob_not
	conn.commit()
time2 = TI.time()	
time_tot = time2 - time_init
# always call this when finishing up
MPI.Finalize()
cur.close()
cur2.close()
cur3.close()
conn.close()
conn2.close()
conn3.close()	
print "Time to do 10 million:", time_tot
