/* Constants for the JPL DE421 Chebyshev Ephemeris
 */

#define DE102F 0
#define DE200F 0
#define DE245F 0
#define DE400F 0
#define DE421BSP 1
#define BSP_FORMAT 1

double JDb = 2414864.5; /* starting JED 1899 Jul 29 */
double JDe = 2471184.5; /* ending JED 2053 Oct 9 */
double JDi = 16.0; /* days per record varies */
long recsiz = 16; /* bytes per record varies */
long rec0 = 4096L; /* offset to first record varies */

struct ephloc
{
  long first_item_ordinal;
  long last_item_ordinal;
  int days_per_record;
  int items_per_record;
  int total_records;
};

struct ephloc objs[15] = {
  {513,    310276, 8, 44, 7040},      /* Mercury */
  {310277, 422920, 16, 32, 3520},  /* Venus */
  {422921, 567244, 16, 41, 3520},  /* EMB */
  {567245, 628848, 32, 35, 1760},  /* Mars */
  {628849, 674612, 32, 26, 1760},
  {674613, 715096, 32, 23, 1760},
  {715097, 750300, 32, 20, 1760},
  {750301, 785504, 32, 20, 1760},
  {785505, 820708, 32, 20, 1760},
  {820709, 943912, 16, 35, 3520},  /* Sun */
  {943913, 1521196, 4, 41, 14080}, /* Moon re EMB */
  {1521197, 2098480, 4, 41, 14080}, /* Earth re EMB */
  {2098481, 2098492, 56320, 8, 1},  /* Mercury point */
  {2098493, 2098504, 56320, 8, 1},  /* Venus point */
  {2098505, 2098516, 56320, 8, 1}  /* Mars point */
#if 0
  {  3, 14, 4},
  {171, 10, 2},
  {231, 13, 2},
  {309, 11, 1},
  {342, 8, 1},
  {366, 7, 1},
  {387, 6, 1},
  {405, 6, 1},
  {423, 6, 1},
  {441, 13, 8},
  {753, 11, 2},
  {819, 10, 4},
  {899, 10, 4}
#endif
};

double au = 1.495978706996262e8; /* Kilometers per au */

double emrat = 81.3005690699153; /* Earth/Moon mass ratio */
double radearth = 6378.1363;
double clight = 2.99792458e5;
char *dename = "de421.bsp";
