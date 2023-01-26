import sympy
import pandas as pd

primes = list(sympy.primerange(30))
prime_gaps = [b - a for a, b in zip(primes, primes[1:])]
print(primes)
print(prime_gaps)
import pandas as pd

df = pd.DataFrame()  # prime_gaps)#, index=True)
df["x"] = list(range(1, len(prime_gaps) + 1))
df["prime_gaps"] = prime_gaps

print(df)

# import statsmodels.api as sm

# logit_model = sm.Logit(df["prime_gaps"], df["x"])
# model_results = logit_model.fit()
# print(model_results)

import numpy

# import scipy
from scipy.optimize import curve_fit
from scipy.stats import f_oneway
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

x_train = numpy.array(list(range(1, len(prime_gaps) + 1)))
y_train = numpy.array(prime_gaps)


model1 = sm.OLS(y_train, x_train)
result = model1.fit()
print(result.summary())


# def model(x, a, b):
#     return a + b * numpy.log(x)


# params, cov = curve_fit(model, x, y)

# print(params)
# print(cov)

# F, p = f_oneway(x, y, params)

# print(F)
# print(p)
