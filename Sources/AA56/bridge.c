

#include "kep.h"
#include "bridge.h"

double respolar[3];
/* approximate motion of right ascension and declination
 * of object, in radians per day
 */
extern double FAR dradt;
extern double FAR ddecdt;

/* Space for star description read from a disc file.
 */
extern struct star fstar;

/* Space for orbit read from a disc file. Entering 99 for the
 * planet number yields a prompt for a file name containg ASCII
 * strings specifying the elements.
 */
extern struct orbit forbit;

/* Orbits for each planet.  The indicated orbital elements are
 * not actually used, since the positions are are now calculated
 * from a formula.  Magnitude and semidiameter are still used.
 */
 /* Programs to compute perturbations. */
extern struct plantbl mer404, ven404, ear404, mar404;
extern struct plantbl jup404, sat404, ura404, nep404, plu404;

extern struct orbit mercury;
extern struct orbit venus;

extern struct orbit earth;
extern struct orbit earth;

extern struct orbit mars;
extern struct orbit jupiter;
extern struct orbit saturn ;

extern struct orbit uranus;

extern struct orbit neptune;

extern struct orbit pluto ;


/* coordinates of object
 */
extern int objnum;	/* I.D. number of object */
extern double robject[3]; /* position */
/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
extern double FAR obpolar[3];

/* coordinates of Earth
 */
/* Heliocentric rectangular equatorial position
 * of the earth at time TDT re equinox J2000
 */
extern double FAR rearth[3];
/* Corresponding polar coordinates of earth:
 * longitude and latitude in radians, radius in au
 */
extern double FAR eapolar[3];

/* Julian date of ephemeris
 */
extern double JD;
extern double TDT;
extern double UT;
extern double deltat_value;

/* flag = 0 if TDT assumed = UT,
 *      = 1 if input time is TDT,
 *      = 2 if input time is UT.
 */
extern int jdflag;

/* correction vector, saved for display  */
extern double dp[3];

/* display formats for printf()
 */
extern char *intfmt, *dblfmt;

/* display enable flag
 */
extern int prtflg;

extern struct orbit *elobject;	/* pointer to orbital elements of object */



struct Polar calcPolar(double UnixTimeStamp, int planet)
{
    struct Polar res;
    int i;


    JD = UnixTimeStamp;
    JD = JD / 86400 + 2440587.5;  // Convert to Julian day
    update(); /* find UT and ET */

    objnum= planet;
    switch(objnum)
        {
        case 0: elobject = 0;
            break;
        case 1: elobject = &mercury; break;
        case 2: elobject = &venus; break;
        case 3: elobject = 0;
            break;
        case 4: elobject = &mars; break;
        case 5: elobject = &jupiter; break;
        case 6: elobject = &saturn; break;
        case 7: elobject = &uranus; break;
        case 8: elobject = &neptune; break;
        case 9: elobject = &pluto; break;
        default:
            return res;
        }

    /* Always calculate heliocentric position of the earth */
    kepler( TDT, &earth, rearth, eapolar );



    switch( objnum )
    {
    case 0:
        dosun();
        break;
    case 3:
        domoon();
        break;
    default:
        doplanet();
    break;
    }

    res.lon = respolar[0];
    res.dec = respolar[1];
    res.r = respolar[2];
    return res;
}


struct Polar calcPolarPath(double UnixTimeStamp, int index, char path[])
{
    struct Polar res;
    strcpy(starnam,path);
    strcpy(orbnam,path);
    int i;

//    kinit();


    JD = UnixTimeStamp;
    JD = JD / 86400 + 2440587.5;  // Convert to Julian day
    update(); /* find UT and ET */

    objnum= index;
    switch(objnum)
        {
        case 88:
            elobject = (struct orbit *)&fstar;
            i = getstar( (struct star *) elobject );
            if( i == 1 )
                return res;
            break;

        case 99:
            elobject = &forbit;
            i = getorbit( elobject );
            if( i == 1)
                return res;
            break;

        default:
            return res;
        }

    /* Always calculate heliocentric position of the earth */
    kepler( TDT, &earth, rearth, eapolar );


    switch( objnum )
    {
    case 88:
        dostar();
        break;
    default:
        doplanet();
    break;
    }

    res.lon = obpolar[0];
    res.dec = obpolar[1];
    res.r = obpolar[2];
    return res;
}
int initCalc(char path[])
{
    return kinit(path);
}

