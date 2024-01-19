import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { DirectionMatchStatusComponent } from './direction-match-status.component';

describe('DirectionMatchStatusComponent', () => {
  let component: DirectionMatchStatusComponent;
  let fixture: ComponentFixture<DirectionMatchStatusComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ DirectionMatchStatusComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(DirectionMatchStatusComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
