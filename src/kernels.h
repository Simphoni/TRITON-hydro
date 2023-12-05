/** @file kernels.h
 *  @brief Header containing the Kernels class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for Kernels class
 *
 *  @author Mario Morales Hernandez
 *  @author Md Bulbul Sharif
 *  @author Tigstu T. Dullo
 *  @author Sudershan Gangrade
 *  @author Alfred Kalyanapu
 *  @author Sheikh Ghafoor
 *  @author Shih-Chieh Kao
 *  @author Katherine J. Evans
 *  @bug No known bugs.
 */



#ifndef KERNELS_H
#define KERNELS_H

#include "constants.h"

namespace Kernels
{

/** @brief It calculates flux in x direction.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param dx Cell size
*  @param dt Time step size
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain (elevation)
*  @param sqrth Square root of water depth array
*  @param rhsh0 RHS array to store partial water depth contributions (east-north)
*  @param rhsh1 RHS array to store partial water depth contributions (west-south)
*  @param rhsqx0 RHS array to store partial discharge in x direction contributions (east-north)
*  @param rhsqx1 RHS array to store partial discharge in x direction contributions (west-south)
*  @param rhsqy0 RHS array to store partial discharge in y direction contributions (east-north)
*  @param rhsqy1 RHS array to store partial discharge in y direction contributions (west-south)
*  @param hextra Minimum depth (tolerance below water is at rest)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void flux_x(int size, int nrows, int ncols, T dx, T dt, T *h_arr, T *qx_arr, T *qy_arr, T *dem, T *sqrth, 
	T *rhsh0, T *rhsh1, T *rhsqx0, T *rhsqx1, T *rhsqy0, T *rhsqy1, T hextra)
	{

		/**************************************
		 *	RHS sketch
		 *			 _____ _____
		 *			|		|     |
		 *			|	  0|1    |
		 *			|_____|_____|
		 *
		 * ************************************/

#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int ix = (id / ncols);	//row id
			int iy = (id % ncols);	//col id
			int rx = (ix * ncols);	//row offset

			bool
			is_top = (ix <= GHOST_CELL_PADDING-1),
			is_btm = (ix >= nrows - GHOST_CELL_PADDING),
			is_lt = (iy <= GHOST_CELL_PADDING-2),
			is_rt = (iy >= ncols - GHOST_CELL_PADDING);

			if (is_top || is_rt || is_lt || is_btm)  //top, right and bottom are not computed. For the X-flux they dont need to be computed. Caution! If domain decomposition changes, this would change (left and right interfaces should be reconsidered)*/
			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}

			int id1, id2;
			T nx, ny;

			id1 = id;
			id2 = rx + iy + 1;
			nx = 1.0;
			ny = 0.0;

			T h1 = h_arr[id1];
			T h2 = h_arr[id2];

			if (h1 > 0.0 || h2 > 0.0)
			{
				int i;
				T alpha[3], beta[3], eigen[3], eigenE[3]={0.0}, eigev[3][3], d1[3]={0.0}, d2[3]={0.0};
				T u1, u2, v1, v2;

				T qx1 = qx_arr[id1];
				T qy1 = qy_arr[id1];
				T z1 = dem[id1];

				T qx2 = qx_arr[id2];
				T qy2 = qy_arr[id2];
				T z2 = dem[id2];

				if (h1 >= hextra)
				{
					u1 = qx1 / h1;
					v1 = qy1 / h1;
				}
				else
				{
					u1 = 0.0;
					v1 = 0.0;
				}
				if (h2 >= hextra)
				{
					u2 = qx2 / h2;
					v2 = qy2 / h2;
				}
				else
				{
					u2 = 0.0;
					v2 = 0.0;
				}

				T sqh1 = sqrth[id1];
				T sqh2 = sqrth[id2];

				T h12 = 0.5*(h1 + h2);
				T dh = h2 - h1;
				T c12 = sqrt(G*h12);

				T u12 = (u1*sqh1 + u2 * sqh2) / (sqh1 + sqh2);
				T v12 = (v1*sqh1 + v2 * sqh2) / (sqh1 + sqh2);

				T un = u12 * nx + v12 * ny;

				eigen[0] = un - c12;
				eigen[1] = un;
				eigen[2] = un + c12;

				//entropy correction

				T e1 = u1 * nx + v1 * ny - SQRTG * sqh1;
				T e2 = u2 * nx + v2 * ny - SQRTG * sqh2;
				if (e1<0.0 && e2>0.0)
				{
					eigenE[0] = eigen[0] - e1 * (e2 - eigen[0]) / (e2 - e1);
					eigen[0] = e1 * (e2 - eigen[0]) / (e2 - e1);
				}

				e1 = u1 * nx + v1 * ny + SQRTG * sqh1;
				e2 = u2 * nx + v2 * ny + SQRTG * sqh2;

				if (e1<0.0 && e2>0.0)
				{
					eigenE[2] = eigen[2] - e2 * (eigen[2] - e1) / (e2 - e1);
					eigen[2] = e2 * (eigen[2] - e1) / (e2 - e1);
				}

				eigev[0][0] = 1.0;
				eigev[0][1] = u12 - c12 * nx;
				eigev[0][2] = v12 - c12 * ny;
				eigev[1][0] = 0.0;
				eigev[1][1] = -c12 * ny;
				eigev[1][2] = c12 * nx;
				eigev[2][0] = 1.0;
				eigev[2][1] = u12 + c12 * nx;
				eigev[2][2] = v12 + c12 * ny;

				alpha[0] = 0.5*(dh - (((qx2 - qx1)*nx + (qy2 - qy1)*ny) - un * dh) / c12);
				alpha[1] = (((qy2 - qy1) - v12 * dh)*nx - ((qx2 - qx1) - u12 * dh)*ny) / c12;
				alpha[2] = 0.5*(dh + (((qx2 - qx1)*nx + (qy2 - qy1)*ny) - un * dh) / c12);

				T dz = z2 - z1;
				T l1 = z1 + h1;
				T l2 = z2 + h2;
				T dzp = dz;
				T hp;

				if (dz >= 0.0)
				{
					hp = h1;
					if (l1 < z2)
					{
						dzp = h1;
					}
				}
				else
				{
					hp = h2;
					if (l2 < z1)
					{
						dzp = -h2;
					}
				}


				beta[0] = 0.5*G*(hp - 0.5*fabs(dzp))*dzp / c12;

				//wet-wet correction
				hp = h1 + alpha[0];
				if (eigen[0] * eigen[2]<0.0 && hp>0.0 && h1 > 0.0 && h2 > 0.0)  //subcritical
				{
					beta[0] = fmax(beta[0], eigen[0]/(eigen[0]+eigenE[2])*(alpha[0] * eigen[0] - h1 * dx*0.5 / dt + eigenE[2]*h2));
					beta[0] = fmin(beta[0], eigen[2]/(eigen[2]+eigenE[0])*(-alpha[2] * eigen[2] + h2 * dx*0.5 / dt + eigenE[0]*h1));
				}

				beta[1] = 0.0;
				beta[2] = -beta[0];

				l1 = hp - beta[0] / eigen[0]; //left intermediate state
				l2 = hp + beta[2] / eigen[2]; //right intermediate state

				if (l2 < -EPS12 && h2 < EPS12)   //send mass contributions to the left if right cell is dry
				{
					d1[0] += (eigen[0] * alpha[0] - beta[0]) + (eigen[2] * alpha[2] - beta[2]);
				}
				else
				{
					if (l1 < -EPS12 && h1 < EPS12)   //send mass contributions to the right if left cell is dry
					{
						d2[0] += (eigen[0] * alpha[0] - beta[0]) + (eigen[2] * alpha[2] - beta[2]);
					}
					else
					{
						for (i = 0; i < 3; i++)   //regular contributions
						{
							if (eigen[i] > 0.0)
							{
								d2[0] += (eigen[i] * alpha[i] - beta[i])*eigev[i][0];
								d2[1] += (eigen[i] * alpha[i] - beta[i])*eigev[i][1];
								d2[2] += (eigen[i] * alpha[i] - beta[i])*eigev[i][2];
							}
							else
							{
								d1[0] += (eigen[i] * alpha[i] - beta[i])*eigev[i][0];
								d1[1] += (eigen[i] * alpha[i] - beta[i])*eigev[i][1];
								d1[2] += (eigen[i] * alpha[i] - beta[i])*eigev[i][2];
							}
						}
					}
				}

				for (i = 0; i < 3; i++)   //entropy
				{
					//entropy
					if (eigenE[i] > 0.0)
					{
						d2[0] += (eigenE[i] * alpha[i])*eigev[i][0];
						d2[1] += (eigenE[i] * alpha[i])*eigev[i][1];
						d2[2] += (eigenE[i] * alpha[i])*eigev[i][2];
					}
					else
					{
						d1[0] += (eigenE[i] * alpha[i])*eigev[i][0];
						d1[1] += (eigenE[i] * alpha[i])*eigev[i][1];
						d1[2] += (eigenE[i] * alpha[i])*eigev[i][2];
					}
				}

				for (i = 0; i < 3; i++)
				{
					if (fabs(d1[i]) < EPS12)
					{
						d1[i] = 0.0;
					}
					if (fabs(d2[i]) < EPS12)
					{
						d2[i] = 0.0;
					}
				}

				rhsh0[id1] = dt * d1[0] / dx;
				rhsqx0[id1] = dt * d1[1] / dx;
				rhsqy0[id1] = dt * d1[2] / dx;

				rhsh1[id2] = dt * d2[0] / dx;
				rhsqx1[id2] = dt * d2[1] / dx;
				rhsqy1[id2] = dt * d2[2] / dx;
			}else{
				rhsh0[id1] = 0.0;
				rhsqx0[id1] = 0.0;
				rhsqy0[id1] = 0.0;

				rhsh1[id2] = 0.0;
				rhsqx1[id2] = 0.0;
				rhsqy1[id2] = 0.0;
			}

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It calculates flux in y direction.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param dx Cell size
*  @param dt Time step size
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain (elevation)
*  @param sqrth Square root of water depth array
*  @param rhsh0 RHS array to store partial water depth contributions (east-north)
*  @param rhsh1 RHS array to store partial water depth contributions (west-south)
*  @param rhsqx0 RHS array to store partial discharge in x direction contributions (east-north)
*  @param rhsqx1 RHS array to store partial discharge in x direction contributions (west-south)
*  @param rhsqy0 RHS array to store partial discharge in y direction contributions (east-north)
*  @param rhsqy1 RHS array to store partial discharge in y direction contributions (west-south)
*  @param hextra Minimum depth (tolerance below water is at rest)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void flux_y(int size, int nrows, int ncols, T dx, T dt, T *h_arr, T *qx_arr, T *qy_arr, T *dem, T *sqrth,
	T *rhsh0, T *rhsh1, T *rhsqx0, T *rhsqx1, T *rhsqy0, T *rhsqy1, T hextra)
	{

		/**************************************
		 *	RHS sketch
		 *			 _____
		 *			|		|
		 *			|	   |
		 *			|__1__|
		 *			|	0	|
		 *			|	   |
		 *			|_____|
		 *
		 * ************************************/

#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int ix = (id / ncols);	//row id
			int iy = (id % ncols);	//col id
			int rx = (ix * ncols);	//row offset

			bool
			is_top = (ix <= GHOST_CELL_PADDING-1),
			is_btm = (ix >= nrows - GHOST_CELL_PADDING+1),
			is_lt = (iy <= GHOST_CELL_PADDING-1),
			is_rt = (iy >= ncols - GHOST_CELL_PADDING);


			if (is_top || is_rt || is_lt || is_btm)   //top, right and left cells are not computed. Bottom is neccesary because we need to compute the N edge, which is equal to the south of rows-1. Caution! If domain decomposition changes, this would change


			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}

			int id1, id2;
			T nx, ny;

			id1 = id;
			id2 = rx - ncols + iy;
			nx = 0.0;
			ny = 1.0;

			T h1 = h_arr[id1];
			T h2 = h_arr[id2];

			if (h1 > 0.0 || h2 > 0.0)
			{
				int i;
				T alpha[3], beta[3], eigen[3], eigenE[3]={0.0}, eigev[3][3], d1[3]={0.0}, d2[3]={0.0};
				T u1, u2, v1, v2;

				T qx1 = qx_arr[id1];
				T qy1 = qy_arr[id1];
				T z1 = dem[id1];

				T qx2 = qx_arr[id2];
				T qy2 = qy_arr[id2];
				T z2 = dem[id2];

				if (h1 >= hextra)
				{
					u1 = qx1 / h1;
					v1 = qy1 / h1;
				}
				else
				{
					u1 = 0.0;
					v1 = 0.0;
				}
				if (h2 >= hextra)
				{
					u2 = qx2 / h2;
					v2 = qy2 / h2;
				}
				else
				{
					u2 = 0.0;
					v2 = 0.0;
				}

				T sqh1 = sqrth[id1];
				T sqh2 = sqrth[id2];

				T h12 = 0.5*(h1 + h2);
				T dh = h2 - h1;
				T c12 = sqrt(G*h12);
				T u12 = (u1*sqh1 + u2 * sqh2) / (sqh1 + sqh2);
				T v12 = (v1*sqh1 + v2 * sqh2) / (sqh1 + sqh2);
				T un = u12 * nx + v12 * ny;

				eigen[0] = un - c12;
				eigen[1] = un;
				eigen[2] = un + c12;

				//entropy correction

				T e1 = u1 * nx + v1 * ny - SQRTG * sqh1;
				T e2 = u2 * nx + v2 * ny - SQRTG * sqh2;
				if (e1<0.0 && e2>0.0)
				{
					eigenE[0] = eigen[0] - e1 * (e2 - eigen[0]) / (e2 - e1);
					eigen[0] = e1 * (e2 - eigen[0]) / (e2 - e1);
				}

				e1 = u1 * nx + v1 * ny + SQRTG * sqh1;
				e2 = u2 * nx + v2 * ny + SQRTG * sqh2;

				if (e1<0.0 && e2>0.0)
				{
					eigenE[2] = eigen[2] - e2 * (eigen[2] - e1) / (e2 - e1);
					eigen[2] = e2 * (eigen[2] - e1) / (e2 - e1);
				}

				eigev[0][0] = 1.0;
				eigev[0][1] = u12 - c12 * nx;
				eigev[0][2] = v12 - c12 * ny;
				eigev[1][0] = 0.0;
				eigev[1][1] = -c12 * ny;
				eigev[1][2] = c12 * nx;
				eigev[2][0] = 1.0;
				eigev[2][1] = u12 + c12 * nx;
				eigev[2][2] = v12 + c12 * ny;

				alpha[0] = 0.5*(dh - (((qx2 - qx1)*nx + (qy2 - qy1)*ny) - un * dh) / c12);
				alpha[1] = (((qy2 - qy1) - v12 * dh)*nx - ((qx2 - qx1) - u12 * dh)*ny) / c12;
				alpha[2] = 0.5*(dh + (((qx2 - qx1)*nx + (qy2 - qy1)*ny) - un * dh) / c12);

				T dz = z2 - z1;
				T l1 = z1 + h1;
				T l2 = z2 + h2;
				T dzp = dz;
				T hp;

				if (dz >= 0.0)
				{
					hp = h1;
					if (l1 < z2)
					{
						dzp = h1;
					}
				}
				else
				{
					hp = h2;
					if (l2 < z1)
					{
						dzp = -h2;
					}
				}

				beta[0] = 0.5*G*(hp - 0.5*fabs(dzp))*dzp / c12;


				//wet-wet correction
				hp = h1 + alpha[0];
				if (eigen[0] * eigen[2]<0.0 && hp>0.0 && h1 > 0.0 && h2 > 0.0)  //subcritical
				{
					beta[0] = fmax(beta[0], eigen[0]/(eigen[0]+eigenE[2])*(alpha[0] * eigen[0] - h1 * dx*0.5 / dt + eigenE[2]*h2));
					beta[0] = fmin(beta[0], eigen[2]/(eigen[2]+eigenE[0])*(-alpha[2] * eigen[2] + h2 * dx*0.5 / dt + eigenE[0]*h1));
				}

				beta[1] = 0.0;
				beta[2] = -beta[0];

				l1 = hp - beta[0] / eigen[0]; //left intermediate state
				l2 = hp + beta[2] / eigen[2]; //right intermediate state

				if (l2 < -EPS12 && h2 < EPS12)   //send mass contributions to the left if right cell is dry
				{
					d1[0] += (eigen[0] * alpha[0] - beta[0]) + (eigen[2] * alpha[2] - beta[2]);
				}
				else
				{
					if (l1 < -EPS12 && h1 < EPS12)   //send mass contributions to the right if left cell is dry
					{
						d2[0] += (eigen[0] * alpha[0] - beta[0]) + (eigen[2] * alpha[2] - beta[2]);
					}
					else
					{
						for (i = 0; i < 3; i++)   //regular contributions
						{
							if (eigen[i] > 0.0)
							{
								d2[0] += (eigen[i] * alpha[i] - beta[i])*eigev[i][0];
								d2[1] += (eigen[i] * alpha[i] - beta[i])*eigev[i][1];
								d2[2] += (eigen[i] * alpha[i] - beta[i])*eigev[i][2];
							}
							else
							{
								d1[0] += (eigen[i] * alpha[i] - beta[i])*eigev[i][0];
								d1[1] += (eigen[i] * alpha[i] - beta[i])*eigev[i][1];
								d1[2] += (eigen[i] * alpha[i] - beta[i])*eigev[i][2];
							}
						}
					}

				}

				for (i = 0; i < 3; i++)   //entropy
				{
					//entropy
					if (eigenE[i] > 0.0)
					{
						d2[0] += (eigenE[i] * alpha[i])*eigev[i][0];
						d2[1] += (eigenE[i] * alpha[i])*eigev[i][1];
						d2[2] += (eigenE[i] * alpha[i])*eigev[i][2];
					}
					else
					{
						d1[0] += (eigenE[i] * alpha[i])*eigev[i][0];
						d1[1] += (eigenE[i] * alpha[i])*eigev[i][1];
						d1[2] += (eigenE[i] * alpha[i])*eigev[i][2];
					}
				}

				for (i = 0; i < 3; i++)
				{
					if (fabs(d1[i]) < EPS12)
					{
						d1[i] = 0.0;
					}
					if (fabs(d2[i]) < EPS12)
					{
						d2[i] = 0.0;
					}
				}

				//in the Y-solver flux we have to accumulate over rhs to NOT overwrite the X-flux values
				rhsh0[id1] += dt * d1[0] / dx;
				rhsqx0[id1] += dt * d1[1] / dx;
				rhsqy0[id1] += dt * d1[2] / dx;

				rhsh1[id2] += dt * d2[0] / dx;
				rhsqx1[id2] += dt * d2[1] / dx;
				rhsqy1[id2] += dt * d2[2] / dx;
			}
#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It calculates each cells new depth, fluxes and updates them.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param dt Time step size
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain (elevation)
*  @param n_arr Manning array in that domain/subdomain
*  @param rhsh0 RHS array to store partial water depth contributions (east-north)
*  @param rhsh1 RHS array to store partial water depth contributions (west-south)
*  @param rhsqx0 RHS array to store partial discharge in x direction contributions (east-north)
*  @param rhsqx1 RHS array to store partial discharge in x direction contributions (west-south)
*  @param rhsqy0 RHS array to store partial discharge in y direction contributions (east-north)
*  @param rhsqy1 RHS array to store partial discharge in y direction contributions (west-south)
*  @param hextra Minimum depth (tolerance below water is at rest)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void update_cells(int size, int nrows, int ncols, T dt, T *h_arr, T *qx_arr, T *qy_arr, T *dem, T *n_arr,
	T *rhsh0, T *rhsh1, T *rhsqx0, T *rhsqx1, T *rhsqy0, T *rhsqy1, T hextra)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int ix = (id / ncols);	//row id
			int iy = (id % ncols);	//col id

			bool
			is_top = (ix <= GHOST_CELL_PADDING-1),
			is_btm = (ix >= nrows - GHOST_CELL_PADDING),
			is_lt = (iy <= GHOST_CELL_PADDING-1),
			is_rt = (iy >= ncols - GHOST_CELL_PADDING);

			if (is_top || is_lt || is_rt || is_btm) //exclude halo cells
			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}
			T hij, qxij, qyij;

			T hn = h_arr[id];

			//update depth
			hij = hn - rhsh0[id] - rhsh1[id];

			//if there are negative depths, there are removed. Should be residual values (close to machine accuracy but negative)
			if (hij < 0.0)
			{
				hij = 0.0;
			}

			//if water is below hextra, velocities are removed
			if (hij < hextra)
			{
				qxij = 0.0;
				qyij = 0.0;
			}
			else
			{
				//update friction and momentum
				T mx = qx_arr[id] - rhsqx0[id] - rhsqx1[id];
				T my = qy_arr[id] - rhsqy0[id] - rhsqy1[id];

				T modM = sqrt(mx*mx + my * my);
				if (n_arr[id] > EPS12 && hn >= hextra && modM > EPS12)
				{
					T tt = dt * G*n_arr[id] * modM / (hn*hn*cbrt(hn));
					qxij = -0.5*(mx - mx * sqrt(1.0 + 4.0*tt)) / tt;
					qyij = -0.5*(my - my * sqrt(1.0 + 4.0*tt)) / tt;
				}
				else
				{
					qxij = mx;
					qyij = my;
				}
			}

			h_arr[id]=hij;
			qx_arr[id]=qxij;
			qy_arr[id]=qyij;

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It updates solid wall condition for all cells except q_y for halo cells.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain (elevation)
*  @param max_h_arr Max water depth array
*  @param hextra Minimum depth (tolerance below water is at rest)
*  @param mpi_tasks Number of MPI tasks (domain decomposition)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void wet_dry(int size, int nrows, int ncols, T dt, T *h_arr, T *qx_arr, T *qy_arr, T *dem, T *max_h_arr, T hextra, int mpi_tasks)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int ix = (id / ncols);	//row id
			int iy = (id % ncols);	//col id
			int rx = (ix * ncols);	//row offset

			bool
			is_top = (ix <= GHOST_CELL_PADDING-1),
			is_btm = (ix >= nrows - GHOST_CELL_PADDING),
			is_lt = (iy <= GHOST_CELL_PADDING-1),
			is_rt = (iy >= ncols - GHOST_CELL_PADDING);


			if (is_top || is_lt || is_rt || is_btm) //exclude halo cells
			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}

			T hij = h_arr[id];
			T zij = dem[id];

			if (hij > hextra)
			{
				if (((hij + zij < dem[rx + iy + 1]) && (h_arr[rx + iy + 1] < EPS12)) || ((hij + zij < dem[rx + iy - 1]) && (h_arr[rx + iy - 1] < EPS12)))
				{
					qx_arr[id] = 0.0;
				}
				if(ix==GHOST_CELL_PADDING && mpi_tasks>1){ //first real row
						if ((hij + zij < dem[rx + ncols + iy]) && (h_arr[rx + ncols + iy] < EPS12))
						{
							qy_arr[id] = 0.0;
						}

				}else{
					if(ix==nrows - 2*GHOST_CELL_PADDING && mpi_tasks>1){ //last real row
						if ((hij + zij < dem[rx - ncols + iy]) && (h_arr[rx - ncols + iy] < EPS12))											  {
							qy_arr[id] = 0.0;
						}
					}else{
						if (((hij + zij < dem[rx - ncols + iy]) && (h_arr[rx - ncols + iy] < EPS12)) || ((hij + zij < dem[rx + ncols + iy]) && (h_arr[rx + ncols + iy] < EPS12)))
						{
							qy_arr[id] = 0.0;
						}
							
					}

				}
			}
			
			if(hij > max_h_arr[id]) max_h_arr[id] = hij;

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It updates q_y for halo cells.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain (elevation)
*  @param max_h_arr Max water depth array
*  @param hextra Minimum depth (tolerance below water is at rest)
*  @param mpi_tasks Number of MPI tasks (domain decomposition)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void wet_dry_qy_halo(int size, int nrows, int ncols, T *h_arr, T *qy_arr, T *dem, T hextra)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif
			int id1, id2;

			id1 =  id+(nrows - 2*GHOST_CELL_PADDING-1)*ncols;

			if(id < ncols*GHOST_CELL_PADDING)
			{
				id1=id;
			}

			id2=id1+ncols;

			T h1 = h_arr[id1];
			T h2 = h_arr[id2];
			T z1 = dem[id1];
			T z2 = dem[id2];

			if (h1 > hextra && h1+z1 < z2 && h2 < EPS12)
			{
				qy_arr[id1] = 0.0;
			}
			if (h2 > hextra && h2+z2 < z1 && h1 < EPS12)
			{
				qy_arr[id2] = 0.0;
			}
			
#ifdef ACTIVE_OMP
		}
#endif
	}




/** @brief It copies all halo cells flux value from halo flux bundle to main discharge array.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param halo_qxqy Bundle array that contains only halo discharge cells
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void halo_copy_to_gpu(int size, int nrows, int ncols, T *h_arr, T *qx_arr, T *qy_arr, T *halo)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int index =  id + 2*ncols*GHOST_CELL_PADDING;
			if(id < ncols*GHOST_CELL_PADDING)
			{
				index = id;
			}
			
			h_arr[id]  = halo[index + 0*GHOST_CELL_PADDING*ncols];
			qx_arr[id] = halo[index + 1*GHOST_CELL_PADDING*ncols];
			qy_arr[id] = halo[index + 2*GHOST_CELL_PADDING*ncols];

			int temp = id + (nrows - 2*GHOST_CELL_PADDING)*ncols;

			h_arr[temp]  = halo[index + 6*GHOST_CELL_PADDING*ncols];
			qx_arr[temp] = halo[index + 7*GHOST_CELL_PADDING*ncols];
			qy_arr[temp] = halo[index + 8*GHOST_CELL_PADDING*ncols];

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It copies all halo cells from main array to halo flux bundle.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param halo_qxqy Bundle array that contains only halo discharge cells
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void halo_copy_from_gpu(int size, int nrows, int ncols, T *h_arr, T *qx_arr, T *qy_arr, T *halo)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int index = id + 2*ncols*GHOST_CELL_PADDING;
			if(id < ncols*GHOST_CELL_PADDING)
			{
				index = id;
			}

			halo[index + 0*GHOST_CELL_PADDING*ncols] = h_arr[id];
			halo[index + 1*GHOST_CELL_PADDING*ncols] = qx_arr[id];
			halo[index + 2*GHOST_CELL_PADDING*ncols] = qy_arr[id];

			int temp = id + (nrows - 2*GHOST_CELL_PADDING)*ncols;

			halo[index + 6*GHOST_CELL_PADDING*ncols] = h_arr[temp];
			halo[index + 7*GHOST_CELL_PADDING*ncols] = qx_arr[temp];
			halo[index + 8*GHOST_CELL_PADDING*ncols] = qy_arr[temp];

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It calculates flow for each flow locations and updates depth value of those cells.
*
*  @param size Array size
*  @param hyd_time Array that contains time information of hydrograph file
*  @param hyd_val Array that contains flow information of hydrograph file
*  @param dx Cell size
*  @param dt Time step size
*  @param simtime Current time of simulation
*  @param _idx_low Time index just before current time
*  @param _idx_high Time index just after current time
*  @param h_arr Water depth array
*  @param pos_arr Flow location position array
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void compute_flow(int size, T *hyd_time, T *hyd_val, T dx, T dt, T simtime, int _idx_low, int _idx_high, T *h_arr, int *pos_arr)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			T flow_at_idx_low = hyd_val[_idx_low*size + id];
			T flow_at_idx_high = hyd_val[_idx_high*size + id];
			T time_at_idx_low = hyd_time[_idx_low];
			T time_at_idx_high = hyd_time[_idx_high];

			T flow;

			if (_idx_low == _idx_high)
			{
				flow = flow_at_idx_low;
			}
			else
			{
				T time_diff = simtime - time_at_idx_low;
				T time_diff_2 = time_at_idx_high - time_at_idx_low;
				flow = flow_at_idx_low + (((flow_at_idx_high - flow_at_idx_low) * time_diff) / time_diff_2);
			}
			int sid = pos_arr[id];
			T hij = h_arr[sid];
			T h_src = (flow * dt) / (dx * dx);
			hij += h_src;
			h_arr[sid]=hij;

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It calculates runoff and updates depth and flux value of those cells.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param dt Time step size
*  @param runoff_id Array that contains runoff id
*  @param index_row_runoff Index of current row runoff
*  @param n_rows_runoff Number of runoff rows
*  @param runoff_intensity Array of runoff intensity
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param hextra Minimum depth (tolerance below water is at rest)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void update_runoff(int size, int nrows, int ncols, T dt, int *runoff_id, int index_row_runoff, int n_rows_runoff, T *runoff_intensity, T *h_arr, T *qx_arr, T *qy_arr, T hextra)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif

			int ix = (id / ncols);	//row id
			int iy = (id % ncols);	//col id

			bool
				is_top = (ix <= GHOST_CELL_PADDING-1),
				is_btm = (ix >= nrows - GHOST_CELL_PADDING),
				is_lt = (iy <= GHOST_CELL_PADDING-1),
				is_rt = (iy >= ncols - GHOST_CELL_PADDING);

			if (is_top || is_lt || is_rt || is_btm) //exclude halo cells
			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}

			int id_r;
			T net;
			T hij = h_arr[id];

			id_r = runoff_id[id];
			//runoff_intensity contains the vector of runoff intensities (nrowsRunoffs*ntimes)
			net = runoff_intensity[id_r*n_rows_runoff + index_row_runoff];
			hij += net * dt;
			//if water is below hextra, velocities are removed
			if (hij < hextra)
			{
				//negative water removed
				if (hij < 0.0)
				{
					hij = 0.0;
				}

				qx_arr[id] = 0.0;
				qy_arr[id] = 0.0;
			}

			h_arr[id]=hij;

#ifdef ACTIVE_OMP
		}
#endif
	}


/** @brief It calculates time step size for each cell. It also performs first step of reduction to find minimum time step size.
*
*  @param size Number of times to execute
*  @param dx Cell size
*  @param input_qx Discharge in x direction array
*  @param input_qy Discharge in y direction array
*  @param input_h Water depth array
*  @param output Time step size holder array
*  @param cn CFL value
*  @param hextra Minimum depth (tolerance below water is at rest)
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void compute_dt_and_sqrt(int size, T dx, T *input_qx, T *input_qy, T *input_h,  T *input_sqrth, T *output, T cn, T hextra)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		int tx = threadIdx.x;
		extern __shared__ T l_output[];

		l_output[tx] = MAX_VALUE;
		if (id < size)
		{
			T hij = input_h[id];
			T sqrthij = fmax(sqrt(hij),0.0);
			input_sqrth[id] = sqrthij;

			if (hij > hextra)
			{				
				T maxu=fmax(fabs(input_qx[id]),fabs(input_qy[id]))/hij;
				l_output[tx]=cn*dx/(maxu+ SQRTG*sqrthij);
			}
		}

		__syncthreads();

		for (int offset = blockDim.x >> 1; offset > 0; offset >>= 1)
		{
			if (tx < offset)
			{
				l_output[tx] = fmin(l_output[tx], l_output[tx + offset]);
			}
			__syncthreads();
		}

		if (tx == 0)
		{
			output[blockIdx.x] = l_output[0];
		}
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
			output[id] = MAX_VALUE;
			T hij = input_h[id];
			T sqrthij = fmax(sqrt(hij),0.0);
			input_sqrth[id] = sqrthij;

			if (hij > hextra)
			{				
				T maxu=fmax(fabs(input_qx[id]),fabs(input_qy[id]))/hij;
				output[id]=cn*dx/(maxu+ SQRTG*sqrthij);
			}
		}
#endif
	}


/** @brief It calculates minimum time step size using further reduction.
*
*  @param size Array size
*  @param input Array that contains all time step sizes
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void find_min_dt(int size, T *input)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		int tx = threadIdx.x;

		extern __shared__ T l_output[];

		if(id >= size)
		{
			l_output[tx] = MAX_VALUE;
		}
		else
		{
			l_output[tx] = input[id];			
		}

		__syncthreads();

		for (int offset = blockDim.x >> 1; offset > 0; offset >>= 1)
		{
			if (tx < offset)
			{
				l_output[tx] = fmin(l_output[tx], l_output[tx + offset]);
			}
			__syncthreads();
		}

		if (tx == 0)
		{
			input[blockIdx.x] = l_output[0];
		}
#else
		T min_val = MAX_VALUE;
#pragma omp parallel for reduction(min : min_val)
		for (int id = 0; id < size; id++)
		{
			min_val = fmin(min_val, input[id]);
		}
		input[0] = min_val;
#endif
	}


/** @brief It calculates external boundary condition in each boundary cells and updates those cells depth.
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param dt Time step size
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param dem DEM array in that domain/subdomain
*  @param n_arr Manning array in that domain/subdomain
*  @param relative_index Index of boundary cells in that subdomain
*  @param extbctype Boundary condition type of each cells
*  @param start_index Start index of each boundary condition
*  @param nrows_vars Number of rows for each variable
*  @param extbcvar1 First variable array of boundary conditions
*  @param extbcvar2 Second variable array of boundary conditions
*  @param simtime Current time of the simulation
*  @param rank Current process number
*  @param total_process Total number of MPI processes
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void compute_extbc_values(int size, int nrows, int ncols, T dt, T *h_arr, T *qx_arr, T *qy_arr, T *dem, T *n_arr, int *relative_index, int *extbctype, int *start_index, int *nrows_vars, T *extbcvar1, T *extbcvar2, T simtime, int rank, int total_process)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
		for (int id = 0; id < size; id++)
		{
#endif
			int ii=relative_index[id];
			int bctype=extbctype[id];
			int start_index2=start_index[id];

			
			int ix = (ii / ncols);	//row id
			int iy = (ii % ncols);	//col id

			bool
			is_top = (ix == 1),
			is_btm = (ix == nrows - 2),
			is_lt = (iy == 1),
			is_rt = (iy == ncols - 2);
			
			T hij,qxij,qyij;

			//get interpolated value
			T auxvalue;
			T lvar;
			int idx_low, idx_high;
			idx_low=-999;
			auxvalue=-999;
			
			if(bctype==1){
				lvar=0.5*(simtime+simtime+dt); //midpoint rule
				idx_high=nrows_vars[id]-1;
				
				if(lvar<extbcvar1[0+start_index2]){
					idx_low = 0;
					idx_high = 0;
				}else{
					if (lvar>extbcvar1[idx_high+start_index2])
					{
						idx_low = idx_high;
					}
					else
					{
						for (int i = 0; i < idx_high; i++)
						{
							if (extbcvar1[i+start_index2] <= lvar && extbcvar1[i+1+start_index2] > lvar )
							{
								idx_low = i;
								idx_high = i + 1;
								break;
							}
						}
					}
				}

				T var1_at_idx_low = extbcvar1[idx_low+start_index2];
				T var1_at_idx_high = extbcvar1[idx_high+start_index2];
				T var2_at_idx_low = extbcvar2[idx_low+start_index2];
				T var2_at_idx_high = extbcvar2[idx_high+start_index2];

				
				if (idx_low == idx_high)
				{
					auxvalue = var2_at_idx_low;
				}
				else
				{
					T time_diff = lvar - var1_at_idx_low;
					T time_diff_2 = var1_at_idx_high - var1_at_idx_low;
					auxvalue = var2_at_idx_low + (((var2_at_idx_high - var2_at_idx_low) * time_diff) / time_diff_2);
				}
			}

			if(bctype==2 || bctype==3){ //normal slope or Froude, only one value--> position 0
				auxvalue=extbcvar1[0+start_index2];
			}

			hij=h_arr[ii];
			qxij=qx_arr[ii];
			qyij=qy_arr[ii];

			if(bctype==1){ //h+z(t)
				hij=fmax(auxvalue-dem[ii],0.0);
			}
			if(bctype==2){ //normal slope
				T vel=sqrt(auxvalue*hij*cbrt(hij)/fmax(n_arr[ii],1e-6));
				qxij=hij*vel;
				qyij=hij*vel;
			}
			if(bctype==3){ //Froude
				T vel=auxvalue*sqrt(G*hij);
				qxij=hij*vel;
				qyij=hij*vel;
			}

			if(is_lt){ //west
				h_arr[ii-1] = hij;
				qx_arr[ii-1] = -qxij;
				qy_arr[ii-1] = 0.0;
				h_arr[ii] = hij;
				qx_arr[ii] = -qxij;
				qy_arr[ii] = 0.0;
			}
			if(is_rt){ //east
				h_arr[ii+1] = hij;
				qx_arr[ii+1] = qxij;
				qy_arr[ii+1] = 0.0;
				h_arr[ii] = hij;
				qx_arr[ii] = qxij;
				qy_arr[ii] = 0.0;
			}
			if (rank == 0 && is_top){ //north
				h_arr[ii-ncols] = hij;
				qx_arr[ii-ncols] = 0.0;
				qy_arr[ii-ncols] = qyij;
				h_arr[ii] = hij;
				qx_arr[ii] = 0.0;
				qy_arr[ii] = qyij;
			}
			if (rank == total_process - 1 && is_btm) //south
			{
				h_arr[ii+ncols] = hij;
				qx_arr[ii+ncols] = 0.0;
				qy_arr[ii+ncols] = -qyij;
				h_arr[ii] = hij;
				qx_arr[ii] = 0.0;
				qy_arr[ii] = -qyij;
			}

#ifdef ACTIVE_OMP
		}
#endif

	}


/** @brief It updates ghost cells for open boundaries
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void copy_info_to_exterior_boundaries_west_east(int size, int nrows, int ncols, T *h_arr, T *qx_arr, T *qy_arr)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif
			
			int id_ghost;
			int id_inner;

			if(id<nrows*GHOST_CELL_PADDING){
				id_ghost=(id%nrows)*ncols+id/nrows;
				id_inner=(id%nrows)*ncols+GHOST_CELL_PADDING;
			}else{
				id_ghost=(id%nrows)*ncols+id/nrows+ncols-(2*GHOST_CELL_PADDING);
				id_inner=(id%nrows)*ncols+ncols-1-GHOST_CELL_PADDING;
			}	
			
			h_arr[id_ghost]=h_arr[id_inner];
			qx_arr[id_ghost]=qx_arr[id_inner];
			qy_arr[id_ghost]=qy_arr[id_inner];

#ifdef ACTIVE_OMP
		}
#endif
	}



/** @brief It updates ghost cells for open boundaries
*
*  @param size Array size
*  @param nrows Number of rows in that domain/subdomain
*  @param ncols Number of columns in that domain/subdomain
*  @param h_arr Water depth array
*  @param qx_arr Discharge in x direction array
*  @param qy_arr Discharge in y direction array
*  @param rank Current process number
*  @param total_process Total number of MPI processes

*/
	template<typename T>
#ifdef ACTIVE_GPU
	__global__
#endif
	void copy_info_to_exterior_boundaries_north_south(int size, int nrows, int ncols, T *h_arr, T *qx_arr, T *qy_arr, int rank, int total_process)
	{
#ifdef ACTIVE_GPU
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if (id >= size)
		return;
#else
#pragma omp parallel for
		for (int id = 0; id < size; id++)
		{
#endif
			
			if (rank > 0 && rank < total_process - 1) //exclude interior domains
			{
#ifdef ACTIVE_GPU
				return;
#else
				continue;
#endif
			}

			int id_ghost;
			int id_inner;

			//it's neccesary to include every assignment inside the if clause to consider the case total_process=1
			if(rank==0){
				id_ghost=id;
				id_inner=(id%ncols)+GHOST_CELL_PADDING*ncols;
				h_arr[id_ghost]=h_arr[id_inner];
				qx_arr[id_ghost]=qx_arr[id_inner];
				qy_arr[id_ghost]=qy_arr[id_inner];
			}

			if(rank == total_process - 1){
				id_ghost=id+(nrows-GHOST_CELL_PADDING)*ncols;
				id_inner=(id%ncols)+(nrows-GHOST_CELL_PADDING-1)*ncols;
				h_arr[id_ghost]=h_arr[id_inner];
				qx_arr[id_ghost]=qx_arr[id_inner];
				qy_arr[id_ghost]=qy_arr[id_inner];
			}


#ifdef ACTIVE_OMP
		}
#endif
	}



}

#endif

