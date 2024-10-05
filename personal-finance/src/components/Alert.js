/* eslint-disable default-case */
import React, { Component } from 'react';
import NumberFormat from 'react-number-format';
import ReactTooltip from 'react-tooltip';

import CONFIG from '../config';

class Alert extends Component {
  constructor(props) {
    super(props);

    let estimated_worth = 0;
    let income_sum = 0;
    let expenses_sum = 0;
    let monthly_revenue = 0;
    let monthly_expense = 0;
    let runway = 0;
    let yearly_deficit = 0;

    props.items.map(function (item) {
      let principal = 0;

      switch (item.type) {
        case CONFIG.ITEM_TYPE.ASSET:
          principal = item.amount;

          estimated_worth += principal;

          monthly_revenue += (item.amount * item.apr) / 100 / 12;
          income_sum += item.amount;
          break;

        case CONFIG.ITEM_TYPE.LIABILITY:
          principal = item.amount;

          estimated_worth -= principal;

          monthly_expense += (item.amount * item.apr) / 100 / 12;
          expenses_sum += item.amount;
          break;

        case CONFIG.ITEM_TYPE.INCOME:
          //principal = ((item.amount * 12) / CONFIG.INCOME_ESTIMATE_APR) * 100;

          //monthly_revenue += item.amount;
          income_sum += item.amount;
          break;

        case CONFIG.ITEM_TYPE.EXPENSE:
          expenses_sum += item.amount;
          break;
      }

      return item;
    });

    if (monthly_revenue < expenses_sum) {
      runway = income_sum / ((expenses_sum - monthly_revenue) * 12);
      yearly_deficit = (expenses_sum - monthly_revenue) * 12;
    }

    this.state = {
      runway: runway,
      estimated_worth: estimated_worth,
      yearly_deficit: yearly_deficit,
    };
  }

  render() {
    return (
      <div className="alert alert-secondary" role="alert">
        <ReactTooltip effect="solid" />
        <p>
          <strong>Net Worth:</strong>&nbsp;
          <span data-tip="Total assets">
            <NumberFormat
              value={this.state.estimated_worth}
              displayType={'text'}
              thousandSeparator={true}
              prefix={'$'}
              decimalScale={0}
            />
          </span>
        </p>
        {/*
        <p>
          <strong>Runway:</strong>&nbsp;
          <span
            data-tip={`Yearly deficit: $${this.state.yearly_deficit.toFixed(
              0
            )}`}
          >
            {this.state.yearly_deficit > 0 && (
              <span>
                <NumberFormat
                  value={this.state.runway}
                  displayType={'text'}
                  thousandSeparator={true}
                  decimalScale={1}
                />{' '}
                Years
              </span>
            )}
            {this.state.yearly_deficit === 0 && <span>&infin;</span>}
          </span>
        </p>
        */}
      </div>
    );
  }
}

export default Alert;
