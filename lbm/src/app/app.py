# Custom imports
from lbm.src.core.factory   import *
from lbm.src.app.huayi import *
from lbm.src.app.poiseuille import *
from lbm.src.app.turek import *
from lbm.src.app.step import *
from lbm.src.app.unit import *
from lbm.src.app.knudsen import *
from lbm.src.app.porous import *
from lbm.src.app.porous_1 import *
from lbm.src.app.porous_kn import *



# Declare factory
app_factory = factory()

# Register apps

app_factory.register("huayi",    huayi)
app_factory.register("poiseuille",    poiseuille)
app_factory.register("turek",    turek)
app_factory.register("step",    step)
app_factory.register("unit",    unit)
app_factory.register("knudsen",    knudsen)
app_factory.register("porous",    porous)
app_factory.register("porous_1",    porous_1)
app_factory.register("porous_kn",    porous_kn)

