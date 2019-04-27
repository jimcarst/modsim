void integrate() {
	double sum_v2 = 0.;

	for (int i = 0; i < n_particles; i++) {
		double r_temp[NDIM] = { 0. };
		for (int d = 0; d < NDIM; ++d) {
			r_temp[d] = 2. * r[i][d] - r_prev_t[i][d] + delta_t * delta_t * f[i][d];
		}
		for (int d = 0; d < NDIM; ++d) {
			r_temp[d] -= floor(r_temp[d] / box[d]) * box[d];//periodic BC
			double delta_r = r_temp[d] - r_prev_t[i][d];
			if (delta_r > 0.5 * box[d]) {// Nearest Image Convention
				delta_r -= box[d];
			}
			if (delta_r < -0.5 * box[d]) {
				delta_r += box[d];
			}
			v[i][d] = delta_r / (2. * delta_t);
			sum_v2 += (v[i][d] * v[i][d]);// / n_particles;
		}
		for (int d = 0; d < NDIM; ++d) {
			r_prev_t[i][d] = r[i][d];
			r[i][d] = r_temp[d];
			//r[i][d] -= floor(r[i][d] / box[d]) * box[d];//periodic BC
		}
	}
	temperature = sum_v2 / (3. * n_particles);
	e_kinetic = 0.5 * sum_v2;
	e_total = (e_potential + 0.5 * sum_v2); //e_total = e_potential + e_kin
	//cout <<"{"<< e_potential << ", " << sum_v2 << "},\n";
}