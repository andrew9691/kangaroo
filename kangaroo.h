#ifndef KANGAROO_H
#define KANGAROO_H

#include <NTL\ZZ.h>
#include <NTL\ZZ_p.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ostream>

NTL::ZZ_p kangaroo(const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& p, const NTL::ZZ& c, const NTL::ZZ& d)
{
	int J = std::floor(std::log2(NTL::conv<double>(d - c)));
	std::vector<NTL::ZZ> S;

	for (std::size_t i = 0; i < J; i++)
		S.emplace_back(NTL::power((NTL::ZZ)2, i));

	size_t n = std::ceil(std::sqrt(NTL::conv<double>(d - c))); // 2*

	std::vector<NTL::ZZ> tame;
	std::vector<NTL::ZZ> wild;

	tame.emplace_back(NTL::PowerMod(a, d, p));
	wild.emplace_back(b);
	std::vector<NTL::ZZ> pT;
	std::vector<NTL::ZZ> pW;
	pT.emplace_back(S[tame[0] % J]);
	pW.emplace_back(S[wild[0] % J]);

	for (std::size_t i = 1; i < n; i++)
	{
		tame.emplace_back(NTL::PowerMod(a, d + pT[i - 1], p));
		pT.emplace_back(pT[i - 1] + S[tame[i] % J]);

		wild.emplace_back((wild[i - 1] * NTL::PowerMod(a, S[wild[i - 1] % J], p)) % p);
		pW.emplace_back(pW[i - 1] + S[wild[i] % J]);
	}

	int times = 0;

m:	NTL::ZZ_p::init(p - 1);
	NTL::ZZ_p res;
	for (std::size_t i = 0; i < tame.size(); i++)
		for (std::size_t j = 0; j < wild.size(); j++)
			if (tame[i] == wild[j])
				return NTL::conv<NTL::ZZ_p>(d + pT[i] - pW[j]);
	
	if (times == 1)
	{
		std::cout << "Not found. Change wild[0]" << std::endl;
		return (NTL::ZZ_p)0;
	}

	int new_n = n + std::floor(n / 2);
	for (std::size_t i = n; i < new_n; i++)
	{
		tame.emplace_back(NTL::PowerMod(a, d + pT[i - 1], p));
		pT.emplace_back(pT[i - 1] + S[tame[i] % J]);

		wild.emplace_back((wild[i - 1] * NTL::PowerMod(a, S[wild[i - 1] % J], p)) % p);
		pW.emplace_back(pW[i - 1] + S[wild[i] % J]);
	}
	times++;
	goto m;
}

#endif
