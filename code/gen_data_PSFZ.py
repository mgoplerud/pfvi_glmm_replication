# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:36:44 2023

@author: Max Goplerud

Load data from P-SF-Z Replication Archive
"""

import os

init_wd = os.getcwd()

os.chdir('crossed_effects_python')

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.special import expit

from xfx.glm.gaussian import sample_posterior as sample_gaussian
from xfx.glm.binomial import sample_posterior as sample_binomial

def generate_PSFZ_data(n_levels_, seed, n_samples, n_warmup, family = 'gaussian'):
    
    def sample_coef_fixture(j, tau, ome):

        alp = [ome.normal(0, 1 / np.sqrt(tau_), j_) for tau_, j_ in zip(tau, j)]
        return [alp_ - np.mean(alp_) for alp_ in alp]


    def sample_randfx_fixture(i, df_tau, scale_tau, ome):
    
        tau = scale_tau * ome.chisquare(df_tau, len(i))
        alp = sample_coef_fixture(i, tau, ome)
        return alp, tau
    
    def sample_mar_design(j, p_miss, ome):
    
        i = np.stack(np.meshgrid(*[np.arange(j_) for j_ in j])).T.reshape(-1, 2)
        i = i[ome.uniform(size=i.shape[0]) > p_miss]
        ome.shuffle(i, 0)
        return i
    
    
    def sample_balanced_fixture(j, alp0=0, df_tau=2, scale_tau=1, ome=np.random.default_rng()):
    
        alp, tau = sample_randfx_fixture(j, df_tau, scale_tau, ome)
        i = sample_balanced_design(j, ome)
        eta = alp0 + np.sum([alp_[j_] for alp_, j_ in zip(alp, i.T)], 0)
        return (eta, i), (alp0, alp, tau)
    
    
    def sample_mar_fixture(j, df_tau=2, scale_tau=1, p_miss=.1, ome=np.random.default_rng()):
    
        alp0 = 0
        alp, tau = sample_randfx_fixture(j, df_tau, scale_tau, ome)
        i = sample_mar_design(j, p_miss, ome)
        eta = alp0 + np.sum([alp_[j_] for alp_, j_ in zip(alp, i.T)], 0)
        return (eta, i), (alp0, alp, tau)

    def package_samples(samples, model, algo, n_obs, run):
    
        bet, tau = zip(*samples)
        alp0 = np.array([bet_[0][0] for bet_ in bet])
        alp = np.array([bet_[1:] for bet_ in bet])
        mean = np.mean(alp, 2).T
        prior_prec = np.array(tau).T
        dfs = [pd.DataFrame({'iter': np.arange(len(samples)), 'value': np.array(alp0), 'factor': [0] * len(samples), 'stat': ['mean'] * len(samples)})]
        for i in range(mean.shape[0]):
            df_mean_ = pd.DataFrame({'iter': np.arange(len(samples)), 'value': mean[i], 'factor': [i + 1] * len(samples), 'stat': ['mean'] * len(samples)})
            df_prior_prec_ = pd.DataFrame({'iter': np.arange(len(samples)), 'value': prior_prec[i], 'factor': [i + 1] * len(samples), 'stat': ['prior_prec'] * len(samples)})
            dfs.extend([df_mean_, df_prior_prec_])
        df = pd.concat(dfs)
        df['model'] = model
        df['algo'] = algo
        df['n_obs'] = n_obs
        df['run'] = run
        return df.set_index(['model', 'algo', 'n_obs', 'run', 'factor', 'stat', 'iter']).unstack('iter').value

    ome = np.random.default_rng(seed)
    fixture = sample_mar_fixture(np.repeat(n_levels_, 2), 1e100, 1e-100, 0.9, ome)[0]
    if family == 'gaussian':        
        gauss_fixture = (ome.normal(fixture[0], 1), 
            None, np.ones_like(fixture[0]), np.repeat(n_levels_, 2), fixture[1])
            
        gauss_collapsed_sampler = sample_gaussian(*gauss_fixture, np.ones(2), np.ones(2), np.inf, 1, None, True, ome)
        gauss_samples = [next(gauss_collapsed_sampler)[:2] for i in range(n_samples + n_warmup)][n_warmup: ]
        out = (gauss_fixture, gauss_samples)
    elif family == 'binomial':
        binom_fixture = (ome.binomial(1, expit(fixture[0])), np.ones_like(fixture[0]),
            np.repeat(n_levels_, 2), fixture[1])
        binom_collapsed_sampler = sample_binomial(*binom_fixture, np.ones(2), np.ones(2), None, True, ome)
        binom_samples = [next(binom_collapsed_sampler)[:2] for i in range(n_samples + n_warmup)][n_warmup: ]
        out = (binom_fixture, binom_samples)
    else:
        raise Exception('Family must be "gaussian" or "binomial"')
    return out

os.chdir(init_wd)
