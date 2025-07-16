### Computational Modeling

### Ordinary Differential Equations (ODE)
Using ODE for Modeling the Imune System on Inflamatory Response. 

Figure shows the diagram for modeling how the tecidual damage behaves in a skin wound without immunologic response. 
![](tp1/modelo1.png)

We can derive Equations as follows: 

$$\frac{d}{dt}TD(t) = \alpha N(t)$$
$$\frac{d}{dt} N(t)= \beta TD(t) + \gamma CH(t) -\alpha N(t)$$
$$\frac{d}{dt} CH(t) = \rho N(t) - \sigma CH(t)$$

Figure includes the immunologic response on the previous model. 
![](tp1/modelo2.png)

$$\frac{d}{dt}TD(t) = \alpha N(t) - u_{reg} M(t)$$
$$\frac{d}{dt}N(t) = \beta TD(t) + \frac{\gamma CH(t)}{(1+\mu_{A}A(t))} -\alpha N(t)$$
$$\frac{d}{dt}CH(t) = \frac{\rho N(t)}{(1+\alpha_A A(t))} -\eta_{CH} CH(t)$$
$$\frac{d}{dt}M(t) = v N(t) -\eta_M M(t)$$
$$\frac{d}{dt}A(t) = w_{reg} M(t) -\eta_A A(t)$$

