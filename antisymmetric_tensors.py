from manim import *
import numpy as np

#FADE OUT EVERYTHING ON SCREEN
def fade_all(self):
    self.play(
        *[FadeOut(mob) for mob in self.mobjects]
    )


class TensorScene(Scene):
    def construct(self):

        Antisymmetric_Tensor = Tex("Anti","symmetric Tensor").scale(1.2)
        self.play(Write(Antisymmetric_Tensor[1].shift(UP*1)))
        self.wait()

        antisymmetric_definition = MathTex("T^{\\mu\\nu} =", "-","T^{\\nu\\mu}")
        construction = MathTex(
            "T",         #0
            "^{\\mu",   #1
            "\\nu}",    #2
            "=",        #3
            "A",        #4
            "^\\mu",    #5
            "B",        #6
            "^\\nu" ,    #7
            "-",        #8
            "A",        #9
            "^\\nu",    #10
            "B",        #11
            "^\\mu"     #12
        )

        self.play(Write(antisymmetric_definition[0]), Write(antisymmetric_definition[2]))
        self.wait()
        self.play(Write(Antisymmetric_Tensor[0].shift(UP*1)), Write(antisymmetric_definition[1]))
        self.wait()
        self.play(
            Antisymmetric_Tensor.animate.shift(UP*1),
            antisymmetric_definition.animate.shift(UP*1)
        )
        self.play(Write(construction))
        self.wait(2)

        inverted = MathTex(
            "T",         #0
            "^{\\nu",   #1
            "\\mu}",    #2
            "=",        #3
            "A",        #4
            "^\\nu",    #5
            "B",        #6
            "^\\mu" ,    #7
            "-",        #8
            "A",        #9
            "^\\mu",    #10
            "B",        #11
            "^\\nu"     #12
        ).shift(DOWN*1)
        construction_copy = construction.copy().shift(DOWN*1)
        self.play(Write(construction_copy))
        self.wait()
        self.play(ReplacementTransform(construction_copy[:3], inverted[:3]))
        self.wait(2)
        self.play(ReplacementTransform(construction_copy[3:], inverted[3:]))

        self.wait(2)
        mu_nu_rect1 = SurroundingRectangle(construction[4:8])
        mu_nu_rect2 = SurroundingRectangle(inverted[9:13])
        mu_nu_rects = VGroup(mu_nu_rect1, mu_nu_rect2)
        nu_mu_rect1 = SurroundingRectangle(construction[9:13], color = GREEN)
        nu_mu_rect2 = SurroundingRectangle(inverted[4:8], color = GREEN)
        nu_mu_rects = VGroup(nu_mu_rect1, nu_mu_rect2)
       
        self.play(Create(mu_nu_rects))
        self.wait()
        self.play(Create(nu_mu_rects))

        self.wait(2)
        self.play(FadeOut(VGroup(construction, inverted, mu_nu_rects, nu_mu_rects)))

        first_derivative = MathTex(
            "{\\partial",   #0
            "T",            #1
            "^{\\mu",       #2
            "\\nu}",        #3
            "\\over",       #4
            "{\\partial",   #5
            "x",            #6
            "^\\mu}}",      #7
            "=",            #8
            "-",            #9
            "{\\partial",   #10
            "T",            #11
            "^{\\nu",       #12
            "\\mu}",        #13
            "\\over",       #14
            "{\\partial",   #15
            "x",            #16
            "^\\mu}}",      #17
        )
        self.wait()
        self.play(
            *[mob.animate.shift(UP*0.2) for mob in self.mobjects]
        )
        self.play(Write(first_derivative))

        em_field_tensor = Matrix(
            [
                [0, "-E_x", "-E_y", "-E_z"],
                ["E_x", 0, "-B_z", "B_y"],
                ["E_y", "B_z", 0, "-B_x"],
                ["E_z", "-B_y", "B_x", 0]
            ]
        )
        self.wait()

        self.play(
            *[mob.animate.shift(RIGHT*2.5) for mob in self.mobjects]
        )
        self.play(Write(em_field_tensor.shift(LEFT*2.5)))
        self.wait()
        self.play(
            *[mob.animate.shift(LEFT*2.5) for mob in self.mobjects],
            FadeOut(em_field_tensor)
        )   

        self.wait()

        second_derivative = MathTex(
            "{\\partial^2", #0
            "T",            #1
            "^{\\mu",       #2
            "\\nu}",        #3
            "\\over",       #4
            "{\\partial",   #5
            "x",            #6
            "^\\mu",        #7
            "\\partial",    #8
            "x^\\nu}}",     #9
            "=",            #10
            "-",            #11
            "{\\partial^2", #12
            "T",            #13
            "^{\\nu",       #14
            "\\mu}",        #15
            "\\over",       #16
            "{\\partial",   #17
            "x",            #18
            "^\\mu",        #19
            "\\partial",    #20
            "x^\\nu}}",     #21
        )

        self.play(Write(second_derivative.shift(DOWN*1.5)))
        self.wait()
        order_irrelevance_rect = SurroundingRectangle(second_derivative[11:22])
        order_irrelevance = Tex("Order of Derivatives Doesn't Matter").shift(DOWN*3).set_color("BLUE")
        self.play(Create(order_irrelevance_rect))
        self.wait()
        self.play(Write(order_irrelevance))
        self.wait()
        nu = MathTex("^\\nu")
        mu = MathTex("^\\mu")
        self.play(
            Swap(
                second_derivative[17:20],
                second_derivative[20:22]
            ),
            ReplacementTransform(second_derivative[14], nu.move_to(second_derivative[15].get_center()+UP*0.03)),
            ReplacementTransform(second_derivative[15], mu.move_to(second_derivative[14].get_center()+DOWN*0.03))
        )
        self.play(
            Swap(
                second_derivative[17:20],
                second_derivative[20:22]  
            )
        )
        self.wait()
        self.wait()
        zero = MathTex("0")
        self.play(FadeOut(order_irrelevance_rect))
        self.play(ReplacementTransform(
                VGroup(second_derivative[11:22], nu, mu),
                zero.move_to(second_derivative[11].get_center()+LEFT*0.1)
            )
        )
        self.wait()

