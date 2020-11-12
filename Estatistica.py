import numpy as np
from scipy.stats import norm, f, chi2
class Stat:
    def __init__(self, A=None, B=None, alfa=0.05):
        self.A = A
        self.B = B
        self.alfa = alfa

    def calcStat(self):
        if self.A:
            self.Na = len(self.A)
            self.Xa = np.mean(self.A)
            self.Va = np.var(self.A)
        if self.B:
            self.Nb = len(self.B)       
            self.Xb = np.mean(self.B)
            self.Vb = np.var(self.B)
        return None

    def setData(self, A=None, B=None):
        if A: self.A = A
        if B: self.B = B
    
    def setData2(self, Xa = None, Xb = None, Va = None, Vb = None, Na = None, Nb=None, plinha = None):
        if Xa: self.Xa = Xa
        if Xb: self.Xb = Xb
        if Va: self.Va = Va
        if Vb: self.Vb = Vb
        if Na: self.Na = Na
        if Nb: self.Nb = Nb
        if plinha: self.pl = plinha

    def setConf(self, alfa):
        self.alfa = alfa

    def getData(self, data='AB'):
        if data == 'A' or data =='a': return self.A
        if data == 'B' or data =='b': return self.B
        if data == 'AB' or data =='ab': return [self.A, self.B]

    def getConf(self):
        return self.alfa

    def IC(self, mode = "difMed"):
        '''Calcula o intervalo de confinaça para um determinado conjunto de 
        dados'''
        self.calcStat()
        if mode == 'difMed':
            Zcrit = norm.ppf(1-self.alfa/2)
            e0 = Zcrit * np.sqrt(self.Va**2 / self.Na + self.Vb**2 / self.Nb)
            dX = self.Xa - self.Xb
            return (dX - e0, dX + e0)
        else:
            print("Erro de entrada")

    def TH(self, mode = "med", cond = None):
        '''Realiza um teste de hipóteses dadas as condições de entrada
            Mode: escolha do tipo de analise que esta sendo feita
            Cond: Condição inicial da comparação (por exemplo: mi=25)
            '''
        self.calcStat()
        if mode == 'med':
            print("ainda nao implementado")

        if mode == "menorVar":
            Quicalc = (self.Na - 1) * self.Va / cond
            Quicrit = chi2.ppf(self.alfa, self.Na-1)
            print("Quicalc = ", Quicalc)
            print("Quicritmin = ", Quicrit)
            if Quicalc < Quicrit:
                print("Rejeito H0")
                return False
            else:
                print("Aceito H0")
                return True

        if mode == "maiorVar":
            Quicalc = (self.Na - 1) * self.Va / cond
            Quicrit = chi2.ppf(1-self.alfa, self.Na-1)
            print("Quicalc = ", Quicalc)
            print("Quicritmax = ", Quicrit)
            if Quicalc < Quicrit:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False

        if mode == "diffVar":
            Fcalc = self.Va/ self.Vb
            Fcritmin = f.ppf(self.alfa/2, self.Na-1, self.Nb-1)
            Fcritmax = f.ppf(1-self.alfa/2, self.Na-1, self.Nb-1)
            print("Fcalc = ", Fcalc)
            print("Fcritmin = ", Fcritmin)
            print("Fcritmax = ", Fcritmax)

            if Fcritmin <= Fcalc <= Fcritmax:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False
    def pValor(self, mode = "media", cond = None):
        '''Calcula o p-valor dadas as condições de entrada
            Mode: escolha do tipo de analise que esta sendo feita
            Cond: Condição inicial da comparação (por exemplo: mi=25)
            '''
        if mode == "media":
            print("Ainda nao implementado")
        
        if mode == "maiorProp":
            # Apenas para tornar o codigo flexivel à entradas com o valor
            # inteiro ao invés de inserir a proporção
            if self.pl>1: self.pl /= self.Na
            if cond > 1: cond /= self.Na

            Zcalc = (self.pl - cond) / np.sqrt(cond * (1 - cond) / self.Na)
            Zcrit = norm.ppf(1-alfa)

            print("TERMINAR")
        else:
            print("Erro de entrada")