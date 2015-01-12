
import psycopg2

conn = psycopg2.connect(host='srv01050.soton.ac.uk', user='frohmaier', password='rates', database='frohmaier')
cur = conn.cursor()


f=open('Dithering_Table_Sorted_2010plus.dat')
fin=f.readline()
fin=f.readline()

for line in f:

	ln=line.split('|')
	ujd=float(ln[0])
	seeing_new=float(ln[1])
	ub1_zp=float(ln[2])
	lmt_mag_new=float(ln[3])
	ccdid=float(ln[4])
	goodpixarea=float(ln[5])
	ra_ul=float(ln[6])
	dec_ul=float(ln[7])
	ra_ur=float(ln[8])
	dec_ur=float(ln[9])
	ra_lr=float(ln[10])
	dec_lr=float(ln[11])
	ra_ll=float(ln[12])
	dec_ll=float(ln[13])
	cur.execute("INSERT INTO obs (ujd, seeing_new, ub1_zp_new, lmt_mag_new, ccdid, good_pixel_area, ra_ul, dec_ul, ra_ur, dec_ur, ra_lr, dec_lr, ra_ll, dec_ll) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",(ujd, seeing_new, ub1_zp, lmt_mag_new, ccdid, goodpixarea, ra_ul, dec_ul, ra_ur, dec_ur, ra_lr, dec_lr, ra_ll, dec_ll))

conn.commit()
cur.close()
conn.close()

print 'Done, check DB'	