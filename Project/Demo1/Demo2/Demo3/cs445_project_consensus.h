/*
 * File: cs445_project.h
 */

#define M_INF -9999.9
#define D 0.1
#define E 0.1
#define R_SIZE 5

#define SEQ_SIZE 500

 typedef struct permutation {
	char* string;
	int number;
} permutation;

typedef struct SA_BWT {
	int length;
	int C_table[5];
	char* bwt_seq;
	int* suffix_array;
	int* BWT_RANK[5];
	double* weighted_seqs[4];
	int pattern_length;
	char* pattern;
	char* rev_pattern;
	int gaps_allowed;
} SA_BWT;

typedef struct S_E {
	int s;
	int e;
} S_E;

 struct node {
	int i;
	char c;
	double* match;
	double* gapX;
	double* gapY;
	double* max;
	struct node * prev;
};
typedef struct alignment_results {
	int s;
	double score;
	char n_seq[500];
	char p_seq[500];
} alignment_results; 

char seq1[] = "dKhGCGTGh";
double seq1_weightA[] = {-0.099971516209232, -1.3238077852806198, -0.6936775326290393, -2.446881481353545, -3.45212725495741, -2.2458000927127615, -3.45212725495741, -3.45212725495741, -0.04065157101526907};
double seq1_weightC[] = {-0.5390426437322652, -0.9408122977662088, 0.30559701304596354, -2.638461867278152, 1.5204265933370762, -2.8904008183704413, -1.99778161750665, -3.4521272549574102, 0.7188463193537253};
double seq1_weightG[] = {0.5515629470419413, 0.07593287521115952, -0.41968749534365474, 1.4902903631711346, -2.566852214827215, 1.5139173985358056, -2.9393315539394993, 1.543926856978898, -0.6541462539811809};
double seq1_weightT[] = {-0.14229496347070403, 0.7478997572708075, 0.398633139408666, -2.0965789199546885, -2.562649476676974, -3.2378417127978043, 1.1819238990670855, -3.45212725495741, -0.45918193656628353};

char seq2[] = "GSCYGRGGRv";
double seq2_weightA[] = {-1.4802685731112213, -2.239022480863647, -2.239022480863647, -0.7819717432508633, -1.0891895105019587, 0.10524098938930977, -1.1258158420815916, -2.239022480863647, -0.37040465955362517, -0.23918906235213086};
double seq2_weightC[] = {-2.239022480863647, 1.0607184552010773, 1.4458505260332277, 0.9424083414554497, -1.3500514510816493, 0.3808845606197496, -1.3500514510816493, -2.239022480863647, -0.4515399316162903, 0.7068471186789703};
double seq2_weightG[] = {1.2427637919032497, 0.1469102767582617, -2.239022480863647, -2.239022480863647, 1.3669000459232796, 0.4905447360147099, 1.327209177073795, 1.4811784923401423, 1.102413401775516, 0.17933051474910025};
double seq2_weightT[] = {-0.4290531133031859, -0.8363243614978533, -1.534909618046164, 0.06069307743471073, -2.239022480863647, -2.239022480863647, -1.4802685731112213, -2.239022480863647, -2.0114782166190244, -1.0802387627226644};

char seq3[] = "RGWWCWbndhGTdCYd";
double seq3_weightA[] = {0.059254189937394085, -1.5061274843997163, 0.353468841199554, 0.9350063726443026, -2.427611349883357, 0.8900071501409318, -1.098468482619457, 0.08749237166895696, 0.29526817199065963, -0.49835197233663536, -2.348113625618381, -0.9838430784125943, -0.16912263318466395, -0.734204707625425, -1.160369451500874, -0.17182379439413833};
double seq3_weightC[] = {-1.5933671294293301, -1.6412584684453462, -3.321705534539866, -2.774027889861662, 1.4902177229870037, -1.5956457510741042, 0.26116640440902356, -0.0727217125913395, -0.7585313146474122, 0.1887920378392305, -1.3998159149907223, -2.269628033754845, -2.162406703882106, 1.3202455954834826, 0.2627048797942628, -0.5306016150345979};
double seq3_weightG[] = {0.9678382330585006, 1.4465069670451478, 0.10755742787165258, -0.7511680898667917, -1.821165355926153, -0.625346154367584, 0.5998783972215237, 0.4457318364127929, 0.3879139817766796, -1.2887427810015533, 1.4631689787181725, -0.733862780847688, 0.3184632367015852, -2.0771798920386773, -1.2072727520608226, -0.08032404996106773};
double seq3_weightT[] = {-1.0702935623025431, -3.321705534539866, 0.1677611555192551, -0.6962086220896072, -2.8655169047529605, -0.7651089490448227, -0.14230802297722217, -0.5835224365654009, -0.35995444695214024, 0.5658890432887309, -2.608517338243837, 0.9695458522465882, 0.41830977831942795, -1.6803919947965895, 0.6792427427244194, 0.40991004738552095};

char seq4[] = "hnnGGWWnddWWGGdbWh";
double seq4_weightA[] = {0.2873931371522786, -0.42432832128709497, -0.6707662784591033, -1.4892943831958096, -2.493269934036155, 0.0796300441207614, 0.26022739230966, 0.20357640732573515, 0.1435223353675144, 0.545750014942865, 0.26022739230966, -0.2267732619381562, -2.493269934036155, -2.493269934036155, -0.2267732619381562, -1.214097358509674, 0.4131446718775913, 0.30070418807420213};
double seq4_weightC[] = {0.21053008629861047, -0.12350529459100401, -0.38049563705128703, -2.4932699340361553, -2.4932699340361553, -1.262542631862601, -0.72726641695214, -0.38049563705128703, -0.9595090257350019, -0.72726641695214, -1.262542631862601, -0.38049563705128703, -1.6994364353432307, -2.4932699340361553, -0.72726641695214, -0.12350529459100401, -1.6994364353432307, 0.23061079481069205};
double seq4_weightG[] = {-0.45657587878506506, -0.38049563705128703, -0.01616579559916477, 1.3833852347583353, 1.305067705970472, -1.6994364353432307, -1.6994364353432307, -0.38049563705128703, -0.38049563705128703, -0.38049563705128703, -0.9595090257350019, -0.9595090257350019, 1.3833852347583353, 1.4790944012175267, 0.16911691231256906, -0.12350529459100401, -0.2437675795650807, -0.5827898045642969};
double seq4_weightT[] = {-0.2726248235913401, 0.5034723324874623, 0.5457500149428635, -1.214097358509674, -0.42432832128709497, 0.6989737868481373, 0.5034723324874623, 0.20357640732573515, 0.41314467187759346, -0.14093275568221372, 0.5034723324874623, 0.6252937251125651, -0.9985203482087784, -1.8702475137568733, 0.36472468754622345, 0.6252937251125651, 0.20357640732573515, -0.24943627006879668};

struct alignment_results results[R_SIZE];




