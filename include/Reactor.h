#ifndef ASTARA_REACTOR_H
#define ASTARA_REACTOR_H

#include <vector>
#include <unique_ptr>



namespace astara{

    class Function;

    class Reactor{

    public:
        Reactor();
        ~Reactor()=default;


    protected:

        // transients
        void dRho();
        void dPower();
        void dFuelTemp();

        // reactivity
        void insertControlRod(double length);
        void injectBoron(double boron_concentration);
        void updateCurrentState();



        /**
         * Most of the time heat transfer coefficient and fuel specific heat can be
         * time and temperature depended. So making it a function with arbitrary number of inputs
         * is a better idea
         */
        std::unique_ptr<Function> _convective_heat_transfer_co_efficient;
        std::unique_ptr<Function> _fuel_specific_heat;

        // temperatures
        double _fuel_temperature;
        double _moderator_temperature;

        // delayed neutron fractions
        const unsigned int _number_of_neutron_groups;
        std::vector<double> _delayed_neutron_const;




    private:

    };
}

#endif //ASTARA_REACTOR_H
