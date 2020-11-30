# REPO Overview: Applying singular value decomposition to different applications including particle positions, image compression, and denoising (image) data.

## Part I: Motion of a particle
An experiment was performed to track the x, y, and z positions of a particle using a sensor from t = 0 to t = 100. The **particle_position.mat** file gives these positions. I used **_singular value decomposition_** to learn more about the motion of the particle and conducted rank-1 and rank-2 approximations of the particle. I then plotted the positions of the particle with the rank-1 approximation and again with the rank-2 approximation.

## Part II: Image compression
I downloaded the llama.jpg image using the MATLAB image processing toolbox and then converted in to a grayscale image before using singular value decomposition to conduct a rank-1, rank-20, and lowest rank approximation to still preserve 90% of the energy of the image.

## Part III: Image denoising
I loaded a noisy image from the listed **NoisyImage.mat** file and used the matrix representing the noisy image to denoise it (gray-scale image). I then performed singular value decomposition of the noisy image and denoised it through using a rank-2 approximation of the image.
