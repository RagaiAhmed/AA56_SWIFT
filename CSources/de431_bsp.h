/* Constants for the JPL DE431 Chebyshev Ephemeris
 */

#define DE102F 0
#define DE200F 0
#define DE245F 0
#define DE400F 0
#define DE430BSP 0
#define DE431BSP 1
#define BSP_FORMAT 1
int defd_part_1 = 0;
int defd_part_2 = 0;
double part_1_JDb = -3100014.5; /* -13199.3415 */
double part_1_JDe = 1721425.5; /* J1.0418 */
double part_2_JDb = 1721424.5; /* 1 AD Jan 02 */
double part_2_JDe = 8000002.5; /* J17190.8487 */

double JDb = -3100014.5; /* -13199.3415 */
double JDe = 1721425.5; /* J1.0418 */
double JDi = 16.0; /* days per record varies */
long recsiz = 16; /* bytes per record varies */
long rec0 = 5120; /* offset to first record in file */

struct ephloc
{
  long first_item_ordinal;
  long last_item_ordinal;
  int days_per_record;
  int items_per_record;
  int total_records;
};

struct ephloc *last_objs = NULL;

/* part 1 */
struct ephloc part_1_objs[14] = {
  {153081635, 179599602,  8, 44, 602681}, /* Mercury */
  {143438719, 153081634, 16, 32, 301341}, /* Venus */
  {131083734, 143438718, 16, 41, 301341}, /* EMB */
  {125810245, 131083733, 32, 35, 150671}, /* Mars */
  {121892795, 125810244, 32, 26, 150671}, /* Jupiter */
  {118427358, 121892794, 32, 23, 150671}, /* Saturn */
  {115413934, 118427357, 32, 20, 150671}, /* Uranus */
  {112400510, 115413933, 32, 20, 150671}, /* Neptune */
  {109387086, 112400509, 32, 20, 150671}, /* Pluto */
  { 98840147, 109387085, 16, 35, 301341}, /* Sun */
  {49420342,   98840146,  4, 41,1205361}, /* Moon re EMB */
  {     537,   49420341,  4, 41,1205361}, /* Earth re EMB */
  {     525,        536, 11100032, 8, 1}, /* Mercury point */
  {     513,        524, 11100032, 8, 1}, /* Venus point */
};

/* part 2 */
struct ephloc part_2_objs[14] = {
  {199345631, 233877846,  8, 44, 784823}, /* Mercury */
  {186788443, 199345630, 16, 32, 392412}, /* Venus */
  {170699547, 186788442, 16, 41, 392412}, /* EMB */
  {163832333, 170699546, 32, 35, 196206}, /* Mars */
  {158730973, 163832332, 32, 26, 196206}, /* Jupiter */
  {154218231, 158730972, 32, 23, 196206}, /* Saturn */
  {150294107, 154218230, 32, 20, 196206}, /* Uranus */
  {146369983, 150294106, 32, 20, 196206}, /* Neptune */
  {142445859, 146369982, 32, 20, 196206}, /* Pluto */
  {128711435, 142445858, 16, 35, 392412}, /* Sun */
  {64355986,  128711434,  4, 41,1569645}, /* Moon re EMB */
  {     537,   64355985,  4, 41,1569645}, /* Earth re EMB */
  {     525,        536, 11100032, 8, 1}, /* Mercury point */
  {     513,        524, 11100032, 8, 1}, /* Venus point */
};

/* Copy data for part 1 or part 2 into this structure.  */
struct ephloc objs[14];

double au = 149597870.700; /* Kilometers per au */
double emrat =  0.813005690741906200e02; /* Earth/Moon mass ratio */
double radearth = 6378.1363;
double clight = 2.99792458e5;
char *dename = "de431_bsp_part-1.bsp";
