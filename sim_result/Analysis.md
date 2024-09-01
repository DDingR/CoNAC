
# 1. 시뮬레이션 세팅
- Friction과 Gravity 존재
	- 그 외 uncertainty와 disturbance는 없음
- 동화형이 만들어준 2-link manipulator model
	- Markus, Elisha D., John T. Agee, and Adisa A. Jimoh. "Trajectory control of a two-link robot manipulator in the presence of gravity and friction." 2013 Africon. IEEE, 2013.
- 전체 토크에 norm ball saturation 있음
	- 넘을 경우 ball에 projection되어 포화적용됨

---

# 2. Control Result

### 2.1. CTRL1 Backstepping
- $k_1=1e2, k_2=1e1$로 선정함
	- backstepping control gain으로 notion 참고
- 외력(마찰과 중력)으로 인하여 성능이 가장 좋지 못함
- 3초 정도에 $\tau_1$이 큰 값을 필요로 함. (타 링크의 질량 등으로 인해)
- 그로 인해 projection된 토크가 들어감으로 성능저하의 원인이 됨

### 2.2. CTRL2 Dixon
- $k_1=1, k_2=1$
- $\alpha = 1e3, \rho=0$로 선정
	- 각각 learning rate, e-modification gain
	- Dixon의 원 논문에서는 e-modification이 없어서 $0$으로 선정함
- 외력을 잘 보상하였음
- 가장 성능이 높음 (초기에 가장 빠르게 수렴)
- 제어입력이 크고(포화되고), 진동함
	- SMC같은 signum 항이나 $-ke$와 같은 항이 없음에도 high-gain control과 유사함
	- 성능이 좋지만 선호되지 못함
- weight norm 또한 크게 커짐

### 2.3. CTRL3 Proposed (no ball cstr)
- $k_1=1, k_2=1$
- $\alpha = 1e3, \beta = [1,1,1,0,0,0,0,0]\cdot 1e-1$
	- $\beta$ 는 penalty term gain
	- $\beta_i=0$은 해당 제약조건이 적용되지 않음을 의미
- 제약조건
	- 각 층의 weight norm: $[30,30,100]$; $\lambda_{1,2,3}$
	- Control Max: (적용 안함); $\lambda_{4,5}$
	- Control Min: (적용 안함); $\lambda_{6,7}$
	- Control norm: (적용 안함); $\lambda_{8}$
- weight norm 제약조건dl 없으면 이론적으로 boundedness증명이 안되어 최소로 필요한 제약조건임
- weight norm 제약조건이 작동하여 norm이 과도히 크지 않음
	- 토크 크기 또한 작아짐
		- 마지막 층 weight matrix의 최대 특이값이 작아지기 때문
		- 이로 인해 초반에 큰 힘을 못내서 성능이 낮음
- 초반 이후에는 제어 성능이 CTRL2와 유사함

### 2.4. CTRL 4 Proposed
- $k_1=1, k_2=1$
- $\alpha = 1e3, \beta = [1,1,1,0,0,0,0,1]\cdot 1e-1$
	- $\beta$ 는 penalty term gain
- 제약조건
	- 각 층의 weight norm: $[30,30,100]$; $\lambda_{1,2,3}$
	- Control Max: (적용 안함); $\lambda_{4,5}$
	- Control Min: (적용 안함); $\lambda_{6,7}$
	- Control norm: $50$; $\lambda_{8}$
- 성능이 CTRL2,3와 유사함
- 하지만 control norm 제약조건으로 크기가 작은 것을 볼 수 있음
- 이로 인해 weight norm 제약조건도 계속 inactive임

---

# 3. 결론
- NN을 이용하여 외력을 적절히 보상함을 확인
- weight norm과 control norm 제약조건으로 비선형 제약조건이 잘 적용되는 것을 보임
- 제약조건으로 인한 성능 저하가 없음을 확인함




