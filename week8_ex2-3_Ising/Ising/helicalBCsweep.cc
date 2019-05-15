void sweep_old() {
	int i;
	int nn, sum, delta;
	for (int k = 0; k < N; k++) {
		// Choose a site
		i = std::uniform_int_distribution<int>(0, N * N)(rng);
		// Calculate the sum of the neighbouring spins
		if ((nn = i + 1) >= N) nn -= N;
		sum = lat.at(nn);
		if ((nn = i - 1) < 0) nn += N;
		sum += lat.at(nn);
		if ((nn = i + N) >= N) nn -= N;
		sum += lat.at(nn);
		if ((nn = i - N) < 0) nn += N;
		sum += lat.at(nn);
		// Calculate the change in energy
		delta = sum * lat.at(i);
		//Decide whether to flip spin
		if (delta <= 0) {
			lat.set_at(i, -lat.at(i));
		}
		else if (random_zero_one(rng) < prob[delta]) {
			lat.set_at(i, -lat.at(i));
		}
	}
}