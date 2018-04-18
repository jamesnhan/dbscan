#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define UNCLASSIFIED -1
#define NOISE -2
#define CORE_POINT 1
#define NOT_CORE_POINT 0

struct point
{
	float x;
	float y;
	int cluster_id;
	int type;
	int num_pts;
} point;

int main(int argc, char* argv[])
{
	unsigned int i, j, k, n, min_pts;
	float max_coord, max_dist;
	struct point* points;

	n = atoi(argv[1]);
	max_coord = atof(argv[2]);
	max_dist = atof(argv[3]);
	min_pts = atoi(argv[4]);
	points = malloc(n * sizeof(point));

	srand(time(NULL));
	for (i = 0; i < n; ++i)
	{
		points[i].x = (int)((float)rand()/(float)(RAND_MAX / max_coord));
		//points[i].y = (float)rand()/(float)(RAND_MAX / max_coord);
		points[i].y = 0.0f;
		points[i].cluster_id = UNCLASSIFIED;
	}

	unsigned int cluster_id = 0;
	for (i = 0; i < n; ++i)
	{
		if (points[i].cluster_id == UNCLASSIFIED)
		{
			unsigned int num_pts = 0;
			unsigned int* neighbors = malloc(n * sizeof(unsigned int));
			for (j = 0; j < n; ++j)
			{
				float dist = sqrt(pow(points[i]. x - points[j].x, 2.0f)
					+ pow(points[i].y - points[j].y, 2.0f));
				if (dist <= max_dist && i != j)
				{
					neighbors[num_pts++] = j;
				}
			}

			points[i].num_pts = num_pts;

			if (num_pts < min_pts)
			{
				points[i].type = NOISE;
			}
			else
			{
				points[i].cluster_id = cluster_id++;
				points[i].type = CORE_POINT;

				for (j = 0; j < num_pts; ++j)
				{
					points[neighbors[j]].cluster_id = points[i].cluster_id;

					unsigned int num_ptsj = 0;
					for (k = 0; k < n; ++k)
					{
						float dist = sqrt(pow(points[neighbors[j]]. x - points[k].x, 2.0f)
							+ pow(points[neighbors[j]].y - points[k].y, 2.0f));
						if (dist <= max_dist && j != k)
						{
							++num_ptsj;
						}
					}

					if (num_ptsj >= min_pts)
					{
						points[neighbors[j]].type = CORE_POINT;
					}
					else
					{
						points[neighbors[j]].type = NOT_CORE_POINT;
					}
					points[neighbors[j]].num_pts = num_ptsj;
				}
			}

			free(neighbors);
		}
	}

	for (i = 0; i < n; ++i)
	{
		printf("(%f, %f): %d ", points[i].x, points[i].y, points[i].cluster_id);
		printf("NUM NEIGHBORS: %d ", points[i].num_pts);
		if (points[i].type == CORE_POINT)
		{
			printf("CORE\n");
		}
		else if (points[i].type == NOT_CORE_POINT)
		{
			printf("NOT CORE\n");
		}
		else
		{
			printf("NOISE\n");
		}
	}

	free(points);
	return 0;
}
