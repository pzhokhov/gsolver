#ifdef _ENABLE_CUDA

#include "fhatha.h"
#include <cuda.h>
//#include <cutil.h>
//#include <cutil_inline.h>

#include "cuda_extmath.h"

__global__ void cuda_fhatha_kernel(float_type* data,  size_t N, float_type* phi_buf, float_type* phi_mult, float_type* jbuf, float_type* x, size_t istride, size_t cstride, size_t dist, size_t pointN);

__global__ void cuda_qdht_kernel(float_type* data, float_type* buf, size_t N, float_type* C, float_type* m1,  size_t istride, size_t cstride, size_t idist, size_t pointN);

void fhatha_runmany_cuda(fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, bool do_extra_transfer, float_type* cuda_buf)
{
  float_type* cuda_data;
  float_type* cuda_phibuf; 
  float_type* cuda_phimult;
  float_type* cuda_j1buf; 
  float_type* cuda_x;
  size_t device_freemem=0, device_totalmem=0;
 
  (cudaMemGetInfo(&device_freemem, &device_totalmem));

  size_t cuda_pointN = exp2(floor(log2((double)device_freemem/p->N/sizeof(f_complex)/3)));
  if (cuda_pointN > nN) cuda_pointN=nN;

  int    blocksize = 64;  if (blocksize > cuda_pointN) blocksize=cuda_pointN;
  size_t cuda_piece_size_ = p->N*cuda_pointN*sizeof(f_complex);

  if (do_extra_transfer) 
  {
   (cudaMalloc((void**)&cuda_data,          cuda_piece_size_));
   (cudaMalloc((void**)&cuda_phibuf,      2*cuda_piece_size_));
  }
  

  (cudaMalloc((void**)&cuda_j1buf,       2*(p->N)*sizeof(f_complex ))); 
  (cudaMalloc((void**)&cuda_phimult,       (p->N)*sizeof(float_type)));
  (cudaMalloc((void**)&cuda_x,             (p->N)*sizeof(float_type)));

  (cudaMemcpy(cuda_j1buf,   p->j1,      2*sizeof(f_complex) *(p->N), cudaMemcpyHostToDevice));
  (cudaMemcpy(cuda_phimult, p->phi_mult,  sizeof(float_type)*(p->N), cudaMemcpyHostToDevice));
  (cudaMemcpy(cuda_x,       p->x,         sizeof(float_type)*(p->N), cudaMemcpyHostToDevice));

  if (do_extra_transfer)
  {
  if (cuda_pointN < nN)
  {
   float_type* data_restrided = (float_type*)malloc_ch(cuda_piece_size_); 
   int newstride = cuda_pointN;
 
   for (long i=0; i<nN; i+= cuda_pointN)
   { 
    for (int ni=0; ni<cuda_pointN; ni++)
    for (int nt=0; nt<p->N; nt++) 
    {
	   data_restrided[ni +  2*nt   *newstride]  = real(data[nt*stride + (i+ni)*dist]);
	   data_restrided[ni + (2*nt+1)*newstride]  = imag(data[nt*stride + (i+ni)*dist]);	

	   //data_restrided[ni%(newstride)+(2*nt   + 2*ni/newstride*p->N)*newstride] = real(data[nt*stride + i*dist]);
	   //data_restrided[ni%(newstride)+(2*nt+1 + 2*ni/newstride*p->N)*newstride] = imag(data[nt*stride + i*dist]);
    }
    (cudaMemcpy(cuda_data, data_restrided, cuda_piece_size_, cudaMemcpyHostToDevice)) ;
    cuda_fhatha_kernel<<< cuda_pointN/blocksize,blocksize >>>(cuda_data, p->N, cuda_phibuf, cuda_phimult, cuda_j1buf, cuda_x, 2*newstride, newstride, 1, cuda_pointN);
    (cudaMemcpy(data_restrided, cuda_data, cuda_piece_size_, cudaMemcpyDeviceToHost));
   
    for (int ni=0; ni<cuda_pointN; ni++)
    for (int nt=0; nt<p->N; nt++) 
    {
	   //data[nt*stride+ni*dist] = f_complex(data_restrided[ni%newstride+(2*nt +   ni/newstride*p->N)*newstride], \
		                                   data_restrided[ni%newstride+(2*nt+1 + ni/newstride*p->N)*newstride]);
	   data[nt*stride+(i+ni)*dist] = f_complex(data_restrided[ni +  2*nt   *newstride],\
		                                       data_restrided[ni + (2*nt+1)*newstride]); 
    }
   }
   free(data_restrided);
  }
  else
  {
	  (cudaMemcpy(cuda_data, data, cuda_piece_size_, cudaMemcpyHostToDevice));
	  cuda_fhatha_kernel<<<cuda_pointN/blocksize, blocksize>>>(cuda_data, p->N, cuda_phibuf, cuda_phimult, cuda_j1buf, cuda_x, 2*stride, 1, 2*dist, cuda_pointN);
	  (cudaMemcpy(data, cuda_data, cuda_piece_size_, cudaMemcpyDeviceToHost));
  }
  

  
  (cudaFree(cuda_phibuf));
  (cudaFree(cuda_data));
  }
  else
  {
        cuda_fhatha_kernel<<<cuda_pointN/blocksize, blocksize>>>((float_type*)data, p->N, cuda_buf, cuda_phimult, cuda_j1buf, cuda_x, 2*stride, 1, 2*dist, cuda_pointN);
  }

  (cudaFree(cuda_j1buf));
  (cudaFree(cuda_phimult));
  (cudaFree(cuda_x));
}


__global__ void cuda_fhatha_kernel(float_type* data, size_t N, float_type* phi_buf, float_type* phi_mult, float_type* jbuf, float_type* x, size_t istride, size_t cstride, size_t idist, size_t pointN)
{
	size_t pi      = blockDim.x*blockIdx.x + threadIdx.x; 
	size_t hstride = pointN; 

	data += pi*idist;
	phi_buf += pi;

	for (size_t nt=0; nt < (N-1); nt++) 
	{
		float_type phi_mult_ = phi_mult[nt]; 
		phi_buf[(2*nt)  *hstride] = phi_mult_*(data[(nt)*istride]          -data[(nt+1)*istride        ]); 
		phi_buf[(2*nt+1)*hstride] = phi_mult_*(data[(nt)*istride+cstride]  -data[(nt+1)*istride+cstride]); 
	}

	phi_buf[(2*N-2)*hstride] = data[(N-1)*istride];
	phi_buf[(2*N-1)*hstride] = data[(N-1)*istride+cstride];

	for (size_t nt=2*N; nt < 4*N; nt++) 
	{
		phi_buf[ nt   *hstride] = 0;
	}

	fft_device_strided(phi_buf, 2*N, -1, hstride);

	for (size_t nt=0; nt<2*N; nt++)
	{
		float_type phi_buf_re = phi_buf[(2*nt  )*hstride];
		float_type phi_buf_im = phi_buf[(2*nt+1)*hstride];
		float_type j1re       = jbuf[2*nt];
		float_type j1im       = jbuf[2*nt+1];
		
		phi_buf[(2*nt)  *hstride] = phi_buf_re*j1re - phi_buf_im*j1im;
		phi_buf[(2*nt+1)*hstride] = phi_buf_im*j1re + phi_buf_re*j1im;
	}
	
	fft_device_strided(phi_buf, 2*N, -1, hstride); 
		
	for (size_t nt=0; nt<N; nt++) 
	{
		data[nt*istride]         = phi_buf[ 2*nt   *hstride]/x[nt];
		data[nt*istride+cstride] = phi_buf[(2*nt+1)*hstride]/x[nt];
	}
}


void qdht_runmany_cuda(qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, bool do_extra_transfer, float_type* cuda_buf)
{
  float_type* cuda_data;
  float_type* cuda_buf_; 
  float_type* cuda_C;
  float_type* cuda_m1; 
  size_t device_freemem=0, device_totalmem=0;
 
  (cudaMemGetInfo(&device_freemem, &device_totalmem));

  size_t cuda_pointN = exp2(floor(log2((double)device_freemem/p->N/sizeof(f_complex)/3)));
  if (cuda_pointN > nN) cuda_pointN=nN;

  int    blocksize = 64;  if (blocksize > cuda_pointN) blocksize=cuda_pointN;
  size_t cuda_piece_size = p->N*cuda_pointN*sizeof(f_complex);


  (cudaMalloc((void**)&cuda_C,        (p->N)*(p->N)*sizeof(float_type ))); 
  (cudaMalloc((void**)&cuda_m1,              (p->N)*sizeof(float_type)));

  (cudaMemcpy(cuda_C,   p->C,      sizeof(float_type)*(p->N)*(p->N), cudaMemcpyHostToDevice));
  (cudaMemcpy(cuda_m1,  p->m1,     sizeof(float_type)*(p->N),        cudaMemcpyHostToDevice));

  if (do_extra_transfer)
  {
  (cudaMalloc((void**)&cuda_data,          cuda_piece_size));
  (cudaMalloc((void**)&cuda_buf_,           cuda_piece_size));

  if (cuda_pointN < nN)
  {
   float_type* data_restrided = (float_type*)malloc_ch(cuda_piece_size); 
   int newstride = cuda_pointN;
 
   for (long i=0; i<nN; i+= cuda_pointN)
   { 
    for (int ni=0; ni<cuda_pointN; ni++)
    for (int nt=0; nt<p->N; nt++) 
    {
	   data_restrided[ni +  2*nt   *newstride]  = real(data[nt*stride + (i+ni)*dist]);
	   data_restrided[ni + (2*nt+1)*newstride]  = imag(data[nt*stride + (i+ni)*dist]);	

	   //data_restrided[ni%(newstride)+(2*nt   + 2*ni/newstride*p->N)*newstride] = real(data[nt*stride + i*dist]);
	   //data_restrided[ni%(newstride)+(2*nt+1 + 2*ni/newstride*p->N)*newstride] = imag(data[nt*stride + i*dist]);
    }
    (cudaMemcpy(cuda_data, data_restrided, cuda_piece_size, cudaMemcpyHostToDevice)) ;
    cuda_qdht_kernel<<< cuda_pointN/blocksize,blocksize >>>(cuda_data, cuda_buf_, p->N,  cuda_C, cuda_m1, 2*newstride, newstride, 1, cuda_pointN);
    (cudaMemcpy(data_restrided, cuda_data, cuda_piece_size, cudaMemcpyDeviceToHost));
   
    for (int ni=0; ni<cuda_pointN; ni++)
    for (int nt=0; nt<p->N; nt++) 
    {
	   //data[nt*stride+ni*dist] = f_complex(data_restrided[ni%newstride+(2*nt +   ni/newstride*p->N)*newstride], \
		                                   data_restrided[ni%newstride+(2*nt+1 + ni/newstride*p->N)*newstride]);
	   data[nt*stride+(i+ni)*dist] = f_complex(data_restrided[ni +  2*nt   *newstride],\
		                                       data_restrided[ni + (2*nt+1)*newstride]); 
    }
   }
   free(data_restrided);
  }
  else
  {
	  (cudaMemcpy(cuda_data, data, cuda_piece_size, cudaMemcpyHostToDevice));
	  cuda_qdht_kernel<<<cuda_pointN/blocksize, blocksize>>>(cuda_data, cuda_buf_, p->N, cuda_C, cuda_m1, 2*stride, 1, 2*dist, cuda_pointN);
	  (cudaMemcpy(data, cuda_data, cuda_piece_size, cudaMemcpyDeviceToHost));
  }

  (cudaFree(cuda_data));
  (cudaFree(cuda_buf_));
  }
  else
  {
        cuda_qdht_kernel<<<cuda_pointN/blocksize, blocksize>>>((float_type*)data, cuda_buf, p->N, cuda_C, cuda_m1, 2*stride, 1, 2*dist, cuda_pointN);
  }

  (cudaFree(cuda_C));
  (cudaFree(cuda_m1)); 
}

__global__ void cuda_qdht_kernel(float_type* data, float_type* buf, size_t N, float_type* C, float_type* m1,  size_t istride, size_t cstride, size_t idist, size_t pointN)
{
	size_t pi      = blockDim.x*blockIdx.x + threadIdx.x; 
	size_t hstride = pointN; 

	data += pi*idist;
  	buf += pi; 

	
	for (size_t i=0; i<N; i++) {data[istride*i]/=m1[i]; data[istride*i+cstride]/=m1[i];}
	for (size_t i=0; i<N; i++) 
	{
  	 buf[2*i*hstride] = 0; buf[(2*i+1)*hstride]=0; 
	 for (size_t j=0; j<N; j++) {float_type cC = C[i+N*j]; buf[2*i*hstride]+= cC*data[istride*j]; buf[(2*i+1)*hstride]+=cC*data[istride*j+cstride];}
	}
	for (size_t i=0; i<N; i++) {data[istride*i]=buf[2*i*hstride]*m1[i]; data[istride*i+cstride]=buf[(2*i+1)*hstride]*m1[i]; }
}



#endif
