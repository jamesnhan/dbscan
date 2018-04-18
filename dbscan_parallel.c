#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define UNCLASSIFIED -1
#define NOISE -2
#define CORE_POINT 1
#define NOT_CORE_POINT 0

struct point
{
	double x;
	double y;
	int cluster_id;
	int type;
	int num_pts;
} point;

double min(double a, double b)
{
	return (a < b) ? a : b;
}

double max(double a, double b)
{
	return (a > b) ? a : b;
}

int main(int argc, char* argv[])
{
	int comm_sz, my_rank;
	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	const int nitems = 5;
	int blocklengths[5] = { 1, 1, 1, 1, 1 };
	MPI_Datatype types[5] = { MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
	MPI_Datatype MPI_POINT;
	MPI_Aint offsets[5];

	offsets[0] = offsetof(struct point, x);
	offsets[1] = offsetof(struct point, y);
	offsets[2] = offsetof(struct point, cluster_id);
	offsets[3] = offsetof(struct point, type);
	offsets[4] = offsetof(struct point, num_pts);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);

	unsigned int i, j, k, l, n, min_pts, rem, sum;
	double max_coord, max_dist;

	n = atoi(argv[1]);
	max_coord = atof(argv[2]);
	max_dist = atof(argv[3]);
	min_pts = atoi(argv[4]);

	int sendcounts[comm_sz];
	int displs[comm_sz];
	struct point data[n];
	struct point points[n];
	struct point centroids[n];
	double radius[n];

	rem = n % comm_sz;
	sum = 0;

	for (i = 0; i < comm_sz; ++i)
	{
		sendcounts[i] = n / comm_sz;
		if (rem > 0)
		{
			sendcounts[i]++;
			rem--;
		}

		displs[i] = sum;
		sum += sendcounts[i];
	}

	if (my_rank == 0)
	{
		srand(time(NULL));
		for (i = 0; i < n; ++i)
		{
			data[i].x = (int)((double)rand()/(double)(RAND_MAX / max_coord));
			data[i].y = 0.0;
			data[i].cluster_id = UNCLASSIFIED;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	MPI_Scatterv(data, sendcounts, displs, MPI_POINT, points, sendcounts[my_rank], MPI_POINT, 0, MPI_COMM_WORLD);

	unsigned int cluster_id = 0;

	for (i = 0; i < sendcounts[my_rank]; ++i)
	{
		if (points[i].cluster_id == UNCLASSIFIED)
		{
			unsigned int num_pts = 0;
			unsigned int neighbors[sendcounts[my_rank]];
			for (j = 0; j < sendcounts[my_rank]; ++j)
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

				centroids[points[i].cluster_id].x = points[i].x;
				centroids[points[i].cluster_id].y = points[i].y;

				for (j = 0; j < num_pts; ++j)
				{
					points[neighbors[j]].cluster_id = points[i].cluster_id;
					centroids[points[i].cluster_id].x += points[neighbors[j]].x;
					centroids[points[i].cluster_id].y += points[neighbors[j]].y;

					unsigned int num_ptsj = 0;
					for (k = 0; k < sendcounts[my_rank]; ++k)
					{
						float dist = sqrt(pow(points[neighbors[j]].x - points[k].x, 2.0)
							+ pow(points[neighbors[j]].y - points[k].y, 2.0));
						if (dist <= max_dist && neighbors[j] != k)
						{
							++num_ptsj;
						}
					}

					if (num_ptsj >= min_pts)
					{
						points[neighbors[j]].type = CORE_POINT;
					}
					else if (num_pts > 0)
					{
						points[neighbors[j]].type = NOT_CORE_POINT;
					}
					points[neighbors[j]].num_pts = num_ptsj;
				}

				centroids[points[i].cluster_id].x /= num_pts + 1;
				centroids[points[i].cluster_id].y /= num_pts + 1;
			}

			radius[points[i].cluster_id] = 0.0;
			for (j = 0; j < num_pts; ++j)
			{
				double r = sqrt(pow(points[neighbors[j]].x - centroids[points[i].cluster_id].x, 2.0)
					+ pow(points[neighbors[j]].y - centroids[points[i].cluster_id].y, 2.0));
				if (r > radius[points[i].cluster_id])
				{
					radius[points[i].cluster_id] = r;
				}
			}
		}
	}

	if (my_rank != 0)
	{
		i = my_rank;
		while (i <= log(comm_sz))
		{
			if (i % 2 == 1)
			{
				MPI_Send(centroids, n, MPI_POINT, i + 1, 0, MPI_COMM_WORLD);
				MPI_Send(&cluster_id, n, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
				MPI_Send(points, n, MPI_POINT, i + 1, 0, MPI_COMM_WORLD);
			}
			else
			{
				struct point other_points[n];
				struct point other_centroids[n];
				unsigned int cluster_count;
				MPI_Recv(other_centroids, n, MPI_POINT, i - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&cluster_count, n, MPI_INT, i - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(other_points, n, MPI_POINT, i - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				for (j = 0; j < cluster_id; ++j)
				{
					
					for (k = 0; k < cluster_count; ++k)
					{
						double d = sqrt(pow(centroids[j].x - other_centroids[k].x, 2.0)
							+ pow(centroids[j].x - other_centroids[k].y, 2.0));
						if (d <= abs(radius[j] - radius[k]))
						{
							// Form one big cluster
							for (l = 0; l < sendcounts[i - 1]; ++l)
							{
								if (other_points[l].cluster_id == k)
								{
									other_points[l].cluster_id = j;
									points[sendcounts[my_rank]].x = other_points[l].x;
									points[sendcounts[my_rank]].y = other_points[l].y;
									points[sendcounts[my_rank]].type = other_points[l].type;
									points[sendcounts[my_rank]].cluster_id = other_points[l].cluster_id;
									sendcounts[my_rank]++;
								}
							}

							centroids[j].x = 0.0;
							centroids[j].y = 0.0;
							for (l = 0; l < sendcounts[my_rank]; ++l)
							{
								if (points[l].cluster_id == j)
								{
									centroids[j].x += points[l].x;
									centroids[j].y += points[l].y;
								}
							}
							centroids[j].x /= sendcounts[my_rank];
							centroids[j].y /= sendcounts[my_rank];
							double r = 0.0;
							for (l = 0; l < sendcounts[my_rank]; ++l)
							{
								if (points[l].cluster_id == j)
								{
									double dist = sqrt(pow(points[l].x - centroids[j].x, 2.0)
										+ pow(points[l].y - centroids[j].y, 2.0));
									if (dist > r)
									{
										r = dist;
									}
								}
							}
							radius[j] = r;
						}
						else if (d > abs(radius[j] - radius[k]) && d < radius[j] + radius[k])
						{
							// Intersect in two points
							struct point mid;
							struct point intersection1;
							struct point intersection2;

							double a = (radius[j] * radius[j] - radius[k] * radius[k] + d * d) / (2.0 * d);
							double h = sqrt(radius[j] * radius[j] - a * a);

							mid.x = centroids[j].x + a * (centroids[k].x - centroids[j].x) / d;
							mid.y = centroids[j].y + a * (centroids[k].y - centroids[j].y) / d;

							intersection1.x = mid.x + h * (centroids[k].y - centroids[j].y) / d;
							intersection1.y = mid.y + h * (centroids[k].x - centroids[j].x) / d;
							intersection2.x = mid.x - h * (centroids[k].y - centroids[j].y) / d;
							intersection2.y = mid.y - h * (centroids[k].x - centroids[j].x) / d;

							// Check all points within the rectangle defined by the two intersections
							int num_pts_considered = 0;
							int pts_considered[n];
							for (l = 0; l < sendcounts[my_rank]; ++l)
							{
								if (points[l].x >= min(intersection1.x, intersection2.x)
									&& points[l].x <= max(intersection1.x, intersection2.x))
								{
									if (points[l].y >= min(intersection1.y, intersection2.y)
										&& points[l].y <= max(intersection1.y, intersection2.y))
									{
										pts_considered[num_pts_considered++] = l;
									}
								}
							}

							int m;
							int num_neighbors[num_pts_considered];
							for (l = 0; l < num_pts_considered; ++l)
							{
								for (m = 0; m < num_pts_considered; ++m)
								{
									if (l != m)
									{
										double dist = sqrt(pow(points[pts_considered[l]].x - points[pts_considered[m]].x, 2.0)
											+ pow(points[pts_considered[l]].y - points[pts_considered[m]].y, 2.0));
										if (dist <= max_dist)
										{
											num_neighbors[l]++;
										}
									}
								}
								if (num_neighbors[l] >= min_pts)
								{
									
								}
							}
						}
						else if (d == radius[j] + radius[k])
						{
							// Touch externally
							struct point mid;

							double a = (radius[j] * radius[j] - radius[k] * radius[k] + d * d) / (2.0 * d);
							double h = sqrt(radius[j] * radius[j] - a * a);

							mid.x = centroids[j].x + a * (centroids[k].x - centroids[j].x) / d;
							mid.y = centroids[j].y + a * (centroids[k].y - centroids[j].y) / d;

							for (l = 0; l < sendcounts[my_rank]; ++l)
							{
								double dist = sqrt(pow(points[l].x - mid.x, 2.0)
									+ pow(points[l].y - mid.y, 2.0));
								if (dist <= max_dist)
								{

								}
							}

							for (l = 0; l < sendcounts[i - 1]; ++l)
							{
								double dist = sqrt(pow(other_points[l].x - mid.x, 2.0)
									+ pow(other_points[l].y - mid.y, 2.0));
								if (dist <= max_dist)
								{

								}
							}
						}
					}
				}

				// Examine noise
			}
			i = pow(2.0, i) - 1;
		}
	
	}

	double end = MPI_Wtime();
	double elapsed = end - start;
	double global_elapsed;
	MPI_Reduce(&elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		printf("Runtime: %f\n", global_elapsed);
	}

	MPI_Type_free(&MPI_POINT);
	MPI_Finalize();
	return 0;
}
