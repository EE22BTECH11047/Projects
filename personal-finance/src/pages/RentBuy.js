import React, { Component } from 'react';
import ReactTooltip from 'react-tooltip';

import Footer from '../components/Footer';

import CONFIG from '../config';

class RentBuy extends Component {
  constructor(props) {
    super(props);

    this.state = {
        rent: '',
        property_price: '',
        loan_interest: '',
        loan_duration: '',
        loan_avans: '',
        monthly_payment: 0,
        total_interest: 0
    }

    this.handleInputChange = this.handleInputChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.type === 'checkbox' ? target.checked : target.value;
    const name = target.name;

    this.setState({
      [name]: value
    });
  }

  handleSubmit(e) {
    e.preventDefault();

    let total_interest = (this.state.property_price - this.state.loan_avans) *
        this.state.loan_interest / 100 * this.state.loan_duration;

    let monthly_payment = ((this.state.property_price - this.state.loan_avans) +
        total_interest) / this.state.loan_duration / 12;

    this.setState({
        monthly_payment: monthly_payment,
        total_interest: total_interest
    })
  }

  render() {
    const { rent, property_price, loan_interest, loan_duration, loan_avans } = this.state;
    const valid = rent.length > 0 && property_price.length > 0;

    return (
      <React.Fragment>
        <ReactTooltip effect="solid" />

        <br />

        <div className="container">
            <form onSubmit={this.handleSubmit}>
                <div className="form-row form-group">
                    <div className="col">
                        <input 
                            className="form-control form-control-lg"
                            type="number"
                            placeholder="Monthly Rent"
                            name="rent"
                            value={this.state.rent}
                            onChange={this.handleInputChange} />
                        <small>eg. $100</small>
                    </div>
                    <div className="col">
                        <input 
                            className="form-control form-control-lg"
                            type="number"
                            placeholder="Property Price"
                            name="property_price"
                            value={this.state.property_price}
                            onChange={this.handleInputChange} />
                        <small>eg. $100000</small>
                    </div>
                </div>
                <div className="form-row form-group">
                    <div className="col">
                        <input 
                            className="form-control form-control-lg"
                            type="text"
                            placeholder="Loan Interest"
                            name="loan_interest"
                            value={this.state.loan_interest}
                            onChange={this.handleInputChange} />
                        <small>eg. 4%</small>
                    </div>
                    <div className="col">
                        <input 
                            className="form-control form-control-lg"
                            type="number"
                            placeholder="Loan Duration"
                            name="loan_duration"
                            value={this.state.loan_duration}
                            onChange={this.handleInputChange} />
                        <small>eg. 30 years</small>
                    </div>
                    <div className="col">
                        <input 
                            className="form-control form-control-lg"
                            type="number"
                            placeholder="Loan Avans"
                            name="loan_avans"
                            value={this.state.loan_avans}
                            onChange={this.handleInputChange} />
                        <small>eg. $1000</small>
                    </div>
                </div>
                <div className="form-row form-group">
                    <div className="col">
                        <button
                            className="btn btn-block btn-primary btn-lg"
                            disabled={!valid}>Calculate</button>
                    </div>
                </div>    
            </form>

            <p>Monthly payment: ${this.state.monthly_payment}</p>
            <p>Total Interest: ${this.state.total_interest}</p>
        </div>

        <Footer />
      </React.Fragment>
    );
  }
}

export default RentBuy;
