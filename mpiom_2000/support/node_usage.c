#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <llapi.h>

char *NodeUsage[] = { "SHARED", "NOT_SHARED", "SLICE_NOT_SHARED" };

int main(int argc, char *argv[])
{
  LL_element *query, *job = NULL, *step = NULL;
  int rc, obj_count, err_code, usage;
  int overwrite = 1;
  char **filter_list;
  char *batch, *jobid, *node_usage;

  batch = getenv("LOADLBATCH");
  if (!batch) {
    fprintf(stderr, "Not in LoadLeveler batch mode.\n");
    exit(-1);
  }

  query = ll_query(JOBS);
  if (!query) {
    fprintf(stderr, "ll_query() failed.\n");
    exit(1);
  }

  jobid = getenv("LOADL_STEP_ID");
  if (!jobid) {
    fprintf(stderr, "LOADL_STEP_ID environment variable not available.\n");
    exit(-1);
  }

  memset(strrchr(jobid, '.'), '\0', 1);

  fprintf(stderr, "%s\n", jobid);

  filter_list = (char **) malloc(2*sizeof(char *));
  filter_list[0] = jobid;
  filter_list[1] = NULL;

  rc = ll_set_request(query, QUERY_JOBID, filter_list, ALL_DATA);
  if (rc) {
    fprintf(stderr, "ll_set_request() return code is non-zero.\n");
    exit(1);
  }

  job = ll_get_objs(query, LL_CM, "", &obj_count, &err_code);
  if (job == NULL) {
    printf("ll_get_objs() returns NULL. Error code = %d\n", err_code);
  }

  rc = ll_get_data(job, LL_JobGetFirstStep, &step);
  if (rc) {
    fprintf(stderr, "ll_get_data() return code is non-zero.\n");
    exit(1);
  }

  rc = ll_get_data(step, LL_StepNodeUsage, &usage);
  if (rc) {
    fprintf(stderr, "ll_get_data() return code is non-zero.\n");
    exit(1);
  }

  node_usage = NodeUsage[usage];

  printf("%s\n", node_usage);

  ll_free_objs(query);
  ll_deallocate(query);

  return (0);
}
