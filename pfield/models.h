typedef struct data_node_struct {
  char *description, *data_description;
  int ndata;
  double *data;
  int nchildren;
  struct data_node_struct *children;
  struct data_node_struct *parent;
} data_node;


typedef struct trj_line_struct {
  int tr_idx, mu_idx, phi_idx;
  double mu, theta, cosphi, phi,angbk, angbmkkmz;
} trj_line;

#define NDATA 4

void initializenode(data_node *nodeptr);
data_node *malloc_data_node(int ndata_node);
int loadtrjfile(char *filename, data_node *nodeptr);
int loadgnufile(char *filename, data_node *nodeptr);
int loadintfile(char *filename, data_node *nodeptr);
void printtree_reset(data_node *nodeptr, int reset);
void printtree(data_node *nodeptr);
void evaltree(data_node *nodeptr, double *args, int nargs, double *res);
