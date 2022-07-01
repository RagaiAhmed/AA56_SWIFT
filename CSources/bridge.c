

#include "kep.h"
#include "bridge.h"


/* approximate motion of right ascension and declination
 * of object, in radians per day
 */
double FAR dradt;
double FAR ddecdt;

/* Space for star description read from a disc file.
 */
struct star fstar;

/* Space for orbit read from a disc file. Entering 99 for the
 * planet number yields a prompt for a file name containg ASCII
 * strings specifying the elements.
 */
struct orbit forbit;

/* Orbits for each planet.  The indicated orbital elements are
 * not actually used, since the positions are are now calculated
 * from a formula.  Magnitude and semidiameter are still used.
 */
 /* Programs to compute perturbations. */
extern struct plantbl mer404, ven404, ear404, mar404;
extern struct plantbl jup404, sat404, ura404, nep404, plu404;

struct orbit mercury;
struct orbit venus;

struct orbit earth;
extern struct orbit earth;

struct orbit mars;
struct orbit jupiter;
struct orbit saturn ;

struct orbit uranus;

struct orbit neptune;

struct orbit pluto ;


/* coordinates of object
 */
int objnum;	/* I.D. number of object */
double robject[3]; /* position */
/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
double FAR obpolar[3];

/* coordinates of Earth
 */
/* Heliocentric rectangular equatorial position
 * of the earth at time TDT re equinox J2000
 */
double FAR rearth[3];
/* Corresponding polar coordinates of earth:
 * longitude and latitude in radians, radius in au
 */
double FAR eapolar[3];

/* Julian date of ephemeris
 */
double JD;
double TDT;
double UT;
extern double deltat_value;

/* flag = 0 if TDT assumed = UT,
 *      = 1 if input time is TDT,
 *      = 2 if input time is UT.
 */
int jdflag;

/* correction vector, saved for display  */
double dp[3];

/* display formats for printf()
 */
extern char *intfmt, *dblfmt;

/* display enable flag
 */
int prtflg;

struct orbit *elobject;	/* pointer to orbital elements of object */



struct Polar calcPolar(double UnixTimeStamp, int planet)
{

    int i;

    kinit();

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
            return;
        }

    /* Always calculate heliocentric position of the earth */
    kepler( TDT, &earth, rearth, eapolar );

    struct Polar res;

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
    strcpy(starnam,path);
    strcpy(orbnam,path);
    int i;

    kinit();


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
                return;
            break;

        case 99:
            elobject = &forbit;
            i = getorbit( elobject );
            if( i == 1)
                return;
            break;

        default:
            return;
        }

    /* Always calculate heliocentric position of the earth */
    kepler( TDT, &earth, rearth, eapolar );

    struct Polar res;

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
